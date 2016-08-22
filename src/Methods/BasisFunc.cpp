#include "BasisFunc.h"
#include <cmath>
#include <iostream>
#include <iomanip>

namespace mpi = boost::mpi;
namespace SSAGES
{

	// Pre-simulation hook.
	void Basis::PreSimulation(Snapshot* snapshot, const CVList& cvs)
	{
        // Open file for writing and allocate derivatives vector
        size_t coeff_size = 1, bin_size = 1;
       
        // For print statements and file I/O, the walker IDs are used
        _mpiid = snapshot->GetWalkerID();

        // Make sure the iteration index is set correctly
        _iteration = 0;

        if(!_grid)
		{
            std::cerr<<"::::::::::::::::::::::::::::::::::::::::::::::::::::::"<<std::endl;
			std::cerr<<"ERROR: Method expected a grid but no grid built."<<std::endl;
            std::cerr<<"Exiting on node ["<<_mpiid<<"]"<<std::endl;
            std::cerr<<"::::::::::::::::::::::::::::::::::::::::::::::::::::::"<<std::endl;
			_world.abort(EXIT_FAILURE);
		}

        // There are a few error messages / checks that are in place with defining CVs and grids
        else
        {
            if(_grid->GetDimension() != cvs.size())
            {
                std::cerr<<"::::::::::::::::::::::::::::::::::::::::::::::::::::::"<<std::endl;
                std::cerr<<"ERROR: Grid dimensions doesn't match number of CVS."<<std::endl;
                std::cerr<<"Exiting on node ["<<_mpiid<<"]"<<std::endl;
                std::cerr<<"::::::::::::::::::::::::::::::::::::::::::::::::::::::"<<std::endl;
                _world.abort(EXIT_FAILURE);
            }
            else if(cvs.size() != _polyords.size())
            {
                std::cout<<cvs.size()<<std::endl;
                std::cout<<_polyords.size()<<std::endl;
                std::cout<<"::::::::::::::::::::::::::::::::::::::::::::::::::::::"<<std::endl;
                std::cout<<"WARNING: The number of polynomial orders is not the same"<<std::endl;
                std::cout<<"as the number of CVs"<<std::endl;
                std::cout<<"The simulation will take the first defined input"<<std::endl;
                std::cout<<"as the same for all CVs. ["<<_polyords[0]<<"]"<<std::endl;
                std::cout<<"::::::::::::::::::::::::::::::::::::::::::::::::::::::"<<std::endl;
               
                //! Resize the polynomial vector so that it doesn't crash here
                _polyords.resize(cvs.size());
                //! And now reinitialize the vector
                for(size_t i = 0; i < cvs.size(); ++i)
                {
                    _polyords[i] = _polyords[0];
                }
            }
        }

        // Setting the number of bins here for simplicity
        _nbins.resize(cvs.size());
        for(size_t i = 0; i < cvs.size(); ++i)
            _nbins[i] = _grid->GetNumPoints()[i];

        // This is to check for non-periodic bounds. It comes into play in the update bias function
        _bounds = true;
                 
        for(size_t i = 0; i < cvs.size(); ++i)
        {
            bin_size   *= _nbins[i];
            coeff_size *= _polyords[i]+1;
        }

		_derivatives.resize(cvs.size());
        _histlocal.resize(bin_size,0);
        _histglobal.resize(bin_size,0);
        _unbias.resize(bin_size,0);
        _coeff_arr.resize(coeff_size,0);

        std::vector<int> idx(cvs.size(), 0);
        std::vector<int> jdx(cvs.size(), 0);
		Map temp_map(idx,0.0);
        
        // Initialize the mapping for the hist function
        for(size_t i = 0; i < bin_size; ++i)
        {
            for(size_t j = 0; j < idx.size(); ++j)
            {
                if(idx[j] > 0 && idx[j] % (_nbins[j]) == 0)
                {
                    if(j != cvs.size() - 1)
                    { 
                        idx[j+1]++;
                        idx[j] = 0;
                    }
                }
                temp_map.map[j] = idx[j];
				temp_map.value  = 0.0; 
            } 
            // Resize histogram vectors to correct size
            _hist.push_back(temp_map);
            idx[0]++;
        }
 
        //Initialize the mapping for the coeff function
        for(size_t i = 0; i < coeff_size; ++i)
        {
            for(size_t j = 0; j < jdx.size(); ++j)
            {
                if(jdx[j] > 0 && jdx[j] % (_polyords[j]+1) == 0)
                {
                    if(j != cvs.size() - 1)
                    { 
                        jdx[j+1]++;
                        jdx[j] = 0;
                    }
                }
                temp_map.map[j] = jdx[j];
				temp_map.value  = 0; 
            }
			_coeff.push_back(temp_map);           
            _coeff[i].value = _coeff_arr[i];
            jdx[0]++;
        }

        //Initialize the look-up table.
        BasisInit(cvs);
	}

	// Post-integration hook.
	void Basis::PostIntegration(Snapshot* snapshot, const CVList& cvs)
	{
        std::vector<double> x(cvs.size(),0);
        std::vector<int> idx(cvs.size(),0);
        int ii = 0;

        /*The binned cv space is updated at every step
         *After a certain number of steps has been passed, the system updates a
         *bias projection based on the visited histogram states
         */
        for(size_t i = 0; i < cvs.size(); ++i)
        {
            x[i] = cvs[i]->GetValue();
        }
       
        if(_bounds)
        {
            // Convert the CV value to its discretized value through the grid tool
            idx = _grid->GetIndices(x);
           
            // Map the grid index to the form of the hist and unbias mapping
            for(size_t i = 0; i < cvs.size(); ++i)
            {
                ii += (idx[i])*std::pow(_nbins[i],i);
            } 

            // The histogram is updated based on the index
            _hist[ii].value++;
    
            // Update the basis projection after a predefined number of steps
            if(snapshot->GetIteration()  % _cyclefreq == 0) {	
                double beta;
                beta = 1.0 / (snapshot->GetTemperature() * snapshot->GetKb());

                // For systems with poorly defined temperature (ie: 1 particle) the user needs to define their own temperature. This is a hack that will be removed in future versions. 

                if(snapshot->GetTemperature() == 0)
                {
                    beta = _temperature;
                    if(_temperature == 0)
                    {
                        std::cout<<std::endl;
                        std::cerr<<"::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::"<<std::endl;
                        std::cerr<<"ERROR: Input temperature needs to be defined for this simulation"<<std::endl;
                        std::cerr<<"Exiting on node ["<<_mpiid<<"]"<<std::endl;
                        std::cout<<"::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::"<<std::endl;
                        exit(EXIT_FAILURE);
                    }
                }
                _iteration+= 1;
                UpdateBias(cvs,beta);
                std::cout<<"Node: ["<<_mpiid<<"]"<<std::setw(10)<<"\tSweep: "<<_iteration<<std::endl;
            }
        }

		// This calculates the bias force based on the existing basis projection.
		CalcBiasForce(cvs);

		// Take each CV and add its biased forces to the atoms using the chain rule
		auto& forces = snapshot->GetForces();
		for(size_t i = 0; i < cvs.size(); ++i)
		{
			auto& grad = cvs[i]->GetGradient();

			/* Update the forces in snapshot by adding in the force bias from each
			 *CV to each atom based on the gradient of the CV.
             */
			for (size_t j = 0; j < forces.size(); ++j) 
				for(size_t k = 0; k < 3; ++k) 
					forces[j][k] += _derivatives[i]*grad[j][k];
		}
	}

	// Post-simulation hook.
	void Basis::PostSimulation(Snapshot*, const CVList&)
	{
	    std::cout<<"Run has finished"<<std::endl;	
	}

    /* The basis set is initialized through the recursive definition.
     *Currently, SSAGES only supports Legendre polyonmials for basis projections
     */
    void Basis::BasisInit(const CVList& cvs)
    {
		for( size_t k = 0; k < cvs.size(); k++)
		{
			size_t ncoeff = _polyords[k]+1;

			std::vector<double> dervs(_nbins[k]*ncoeff,0);
			std::vector<double> vals(_nbins[k]*ncoeff,0);
            std::vector<double> x(_nbins[k],0);

            /*As the values for Legendre polynomials can be defined recursively, \
             *both the derivatives and values are defined at the same time,
             */
			for (size_t i = 0; i < _nbins[k]; ++i)
			{
                x[i] = 2.0*(i + 0.5)/(double)(_nbins[k]) - 1.0;
				vals[i] = 1.0;
				dervs[i] = 0.0;
			}

			for (size_t i = 0; i < _nbins[k]; ++i)
			{
				vals[i+_nbins[k]] = x[i];
				dervs[i+_nbins[k]] = 1.0;
			}

			for (size_t j = 2; j < ncoeff; j++)
			{
				for (size_t i = 0; i < _nbins[k]; i++)
				{
                    //Evaluate the values of the Legendre polynomial at each bin
					vals[i+j*_nbins[k]] = ((double) ( 2*j - 1 ) * x[i] * vals[i+(j-1)*_nbins[k]]
					- (double) (j - 1) * vals[i+(j-2)*_nbins[k]]) / (double) (j);

                    //Evaluate the derivatives of the Legendre polynomial at each bin
                    dervs[i+j*_nbins[k]] = ((double) ( 2*j - 1 ) * ( vals[i+(j-1)*_nbins[k]] + x[i] * dervs[i+(j-1)*_nbins[k]] )
                    - (double) (j - 1) * dervs[i+(j-2)*_nbins[k]]) / (double) (j);
				}
			}
            BasisLUT TempLUT(vals,dervs);
            _LUT.push_back(TempLUT);
        }
	}
    
	// Update the coefficients/bias projection
	void Basis::UpdateBias(const CVList& cvs, const double beta)
	{
        std::vector<double> x(cvs.size(), 0);
        std::vector<double> coeffTemp(_coeff.size(), 0);
        double sum  = 0.0;
        double bias = 0.0;
        double basis = 1.0;

        // For multiple walkers, the struct is unpacked
        for(size_t i = 0; i < _hist.size(); ++i)
            _histlocal[i] = (int)_hist[i].value;

        // Summed between all walkers
        MPI_Allreduce(&_histlocal[0], &_histglobal[0], _hist.size(), MPI_INT, MPI_SUM, _world);

        // And then it is repacked into the struct
        for(size_t i = 0; i < _hist.size(); ++i)
            _hist[i].value = _histglobal[i];

        // Construct the biased histogram
        for(size_t i = 0; i < _hist.size(); ++i)
        {
            auto& hist = _hist[i];

            // This is to make sure that the CV projects across the entire surface
            if(hist.value == 0) {hist.value = 1;} 
           
            // The loop builds the previous basis projection for each bin of the histogram
            for(size_t k = 1; k < _coeff.size(); ++k)
            {
                auto& coeff = _coeff[k];
                for(size_t l = 0; l < cvs.size(); ++l)
                { 
                    // The previous bias is only calculated after each sweep has happened
                    basis *= _LUT[l].values[hist.map[l] + coeff.map[l]*(_nbins[l])];
                }
                bias += coeff.value*basis;
                basis = 1.0;
            }
            
            /* The evaluation of the biased histogram which projects the histogram to the
             * current bias of CV space.
             */
            _unbias[i] += hist.value * exp(bias) * _weight / (double)(_cyclefreq); 
            bias = 0.0;
        }

        // The coefficients and histograms are reset after evaluating the biased histogram values
        for(size_t i = 0; i < _coeff.size(); ++i)
        {
            coeffTemp[i] = _coeff[i].value;
            _coeff[i].value = 0.0;
        }

        for(size_t i = 0; i < _hist.size(); ++i)
        {
            _hist[i].value = 0.0;
            _histlocal[i] = 0;
            _histglobal[i] = 0;
        }

        // The loop that evaluates the new coefficients by integrating the CV space
        for(size_t i = 1; i < _coeff.size(); ++i)
        {
            auto& coeff = _coeff[i];
            
            // The method uses a standard integration with trap rule weights
            for(size_t j = 0; j < _hist.size(); ++j)
            {
                auto& hist = _hist[j];
                double weight = std::pow(2.0,cvs.size());

                // This adds in a trap-rule type weighting which lowers error significantly at the boundaries
                for(size_t k = 0; k < cvs.size(); ++k)
                {
                    if(hist.map[k] == 0 || hist.map[k] == _nbins[k]-1)
                        weight /= 2.0;
                }
            
                /*The numerical integration of the biased histogram across the entirety of CV space
                 *All calculations include the normalization as well
                 */
                for(size_t l = 0; l < cvs.size(); l++)
                {
                    basis *= _LUT[l].values[hist.map[l] + coeff.map[l]*(_nbins[l])];
                    basis *=  1.0 / (_nbins[l])*(2 * coeff.map[l] + 1.0);
                }
                coeff.value += basis * log(_unbias[j]) * weight/std::pow(2.0,cvs.size());
                basis = 1.0;
            }
            coeffTemp[i] -= coeff.value;
            _coeff_arr[i] = coeff.value;
            sum += coeffTemp[i]*coeffTemp[i];
        }

        if(_world.rank() == 0)
            // Write coeff at this step, but only one walker
            PrintBias(cvs,beta);

        // The convergence tolerance and whether the user wants to exit are incorporated here
        if(sum < _tol)
        {
            std::cout<<"System has converged"<<std::endl;
            if(_converge_exit)
            {
                std::cout<<"User has elected to exit. System is now exiting"<<std::endl;
                exit(EXIT_SUCCESS);
            }
        }
	}

    /*The coefficients are printed out for the purpose of saving the free energy space
     *Additionally, the current basis projection is printed so that the user can view
     *the current free energy space
     */
    void Basis::PrintBias(const CVList& cvs, const double beta)
    {
        std::vector<double> bias(_hist.size(), 0);
        std::vector<double> x(cvs.size(), 0);
        double temp = 1.0; 
        double pos = 0;

        /* Since the coefficients are the only piece that needs to be
         *updated, the bias is only evaluated when printing
         */
        for(size_t i = 0; i < _hist.size(); ++i)
        {
            for(size_t j = 1; j < _coeff.size(); ++j)
            {
                for(size_t k = 0; k < cvs.size(); ++k)
                {
                    
                    temp *=  _LUT[k].values[_hist[i].map[k] + _coeff[j].map[k] * (_nbins[k])];
                }
                bias[i] += _coeff[j].value*temp;
                temp  = 1.0;
            }
        }

        // The filenames will have a standard name, with a user-defined suffix
        std::string filename1 = "basis"+_bnme+".out";
        std::string filename2 = "coeff"+_cnme+".out";
    
		_basisout.precision(5);
        _coeffout.precision(5);
        _basisout.open(filename1.c_str());
        _coeffout.open(filename2.c_str());

        // The CV values, PMF projection, PMF, and biased histogram are output for the user
        _coeffout << _iteration  <<std::endl;
        _basisout << "CV Values" << std::setw(35*cvs.size()) << "Basis Set Bias" << std::setw(35) << "PMF Estimate" << std::setw(35) << "Biased Histogram" << std::endl;
        
        for(size_t j = 0; j < _unbias.size(); ++j)
        {
            for(size_t k = 0; k < cvs.size(); ++k)
            {
                // Evaluate the CV values for printing purposes
                pos = (_hist[j].map[k]+0.5)*(_grid->GetUpper()[k] - _grid->GetLower()[k]) * 1.0 /(double)( _nbins[k]) + _grid->GetLower()[k];
                _basisout << pos << std::setw(35);
            }
            _basisout << -bias[j] << std::setw(35);
            if(_unbias[j])
                _basisout << -log(_unbias[j]) / beta << std::setw(35);
            else
                _basisout << "0" << std::setw(35);
            _basisout << _unbias[j];
            _basisout << std::endl;
        }

        for(size_t k = 0; k < _coeff.size(); ++k)
        {
            _coeffout <<_coeff[k].value << std::endl;
        }

		_basisout << std::endl;
        _basisout.close();
        _coeffout.close();
	}

    // The forces are calculated by chain rule, first  the derivatives of the basis set, then in the PostIntegration function, the derivative of the CV is evaluated
	void Basis::CalcBiasForce(const CVList& cvs)
	{	
		// Reset derivatives
        std::fill(_derivatives.begin(), _derivatives.end(), 0);
        std::vector<double> x(cvs.size(),0);
        std::vector<int> idx(cvs.size(),0);

        double temp = 1.0;
        size_t ii = 0;

        //This is calculating the derivatives for the bias force
        for (size_t j = 0; j < cvs.size(); ++j)
        {
            x[j] = cvs[j]->GetValue();
            double min = _grid->GetLower()[j];
            double max = _grid->GetUpper()[j];

            if(!_grid->GetPeriodic()[j])
            {
                // In order to prevent the index for the histogram from going out of bounds a check is in place
                if(x[j] > max && _bounds)
                {
                    std::cout<<"::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::"<<std::endl;
                    std::cout<<"WARNING: CV is above the maximum boundary."<<std::endl;
                    std::cout<<"Statistics will not be gathered during this interval"<<std::endl;
                    std::cout<<"::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::"<<std::endl;
                    _bounds = false;
                }
                else if(x[j] < min && _bounds)
                {
                    std::cout<<"::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::"<<std::endl;
                    std::cout<<"WARNING: CV is below the minimum boundary."<<std::endl;
                    std::cout<<"Statistics will not be gathered during this interval"<<std::endl;
                    std::cout<<"::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::"<<std::endl;
                    _bounds = false;
                }
                else if(x[j] < max && x[j] > min && !_bounds)
                {
                    std::cout<<"::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::"<<std::endl;
                    std::cout<<"CV has returned in between bounds. Run is resuming"<<std::endl;
                    std::cout<<"::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::"<<std::endl;
                    _bounds = true;
                }
            }
        }

        // Only apply soft wall potential in the event that it has left the boundaries
        if(_bounds)
        {
            idx = _grid->GetIndices(x);

            for(size_t i = 0; i < cvs.size(); ++i)
                ii += (idx[i])*std::pow(_nbins[i],i);
            
            for (size_t i = 1; i < _coeff.size(); ++i)
            {
                for (size_t j = 0; j < cvs.size(); ++j)
                {
                    temp = 1.0;
                    for (size_t k = 0; k < cvs.size(); ++k)
                    {
                        temp *= j == k ?  _LUT[k].derivs[_hist[ii].map[k] + _coeff[i].map[k]*(_nbins[k])] * 2.0 / (_grid->GetUpper()[j] - _grid->GetLower()[j])
                                       :  _LUT[k].values[_hist[ii].map[k] + _coeff[i].map[k]*(_nbins[k])];
                    }
                    _derivatives[j] -= _coeff[i].value * temp;
                }
            }
        }
        
        // This is where the wall potentials are going to be thrown into the method if the system is not a periodic CV
        for(size_t j = 0; j < cvs.size(); ++j)
        {
            double min = _grid->GetLower()[j];
            double max = _grid->GetUpper()[j];

            if(!_grid->GetPeriodic()[j]) 
            {
                if(x[j] > _boundUp[j])
                    _derivatives[j] -= _restraint[j] * (x[j] - _boundUp[j]);
                else if(x[j] < _boundLow[j])
                    _derivatives[j] -= _restraint[j] * (x[j] - _boundLow[j]);
            }
        }
    }
}

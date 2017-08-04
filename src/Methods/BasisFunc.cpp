/**
 * This file is part of
 * SSAGES - Suite for Advanced Generalized Ensemble Simulations
 *
 * Copyright 2016 Joshua Moller <jmoller@uchicago.edu>
 *           2017 Julian Helfferich <julian.helfferich@gmail.com>
 *
 * SSAGES is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * SSAGES is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with SSAGES.  If not, see <http://www.gnu.org/licenses/>.
 */
#include "BasisFunc.h"
#include "Validator/ObjectRequirement.h"
#include "Drivers/DriverException.h"
#include "CVs/CVManager.h"
#include "Snapshot.h"
#include <cmath>
#include <iostream>
#include <iomanip>

using namespace Json; 

namespace SSAGES
{

	// Pre-simulation hook.
	void BFS::PreSimulation(Snapshot* snapshot, const CVManager& cvmanager)
	{
        auto cvs = cvmanager.GetCVs(cvmask_);
        // For print statements and file I/O, the walker IDs are used
        mpiid_ = snapshot->GetWalkerID();

        // Make sure the iteration index is set correctly
        iteration_ = 0;

        // There are a few error messages / checks that are in place with
        // defining CVs and grids
        if(h_->GetDimension() != cvs.size())
        {
            std::cerr<<"ERROR: Histogram dimensions doesn't match number of CVS."<<std::endl;
            MPI_Abort(world_, EXIT_FAILURE);
        }

        // This is to check for non-periodic bounds. It comes into play in the update bias function
        bounds_ = true;        
        unbias_.resize(h_->size(),0);
 
        // Initialize the mapping for the hist function
        for (auto &val : *h_) {
            val = 0;
        }

        // Reset the values in the force grid
        for (auto &val : *f_) {
            for (size_t i = 0; i <val.size(); i++) {
                val[i] = 0;
            }
        }

        // Reset the values in the bias potential grid
        for (auto &val : *b_) {
            val = 0;
        }
	}

	// Post-integration hook.
	void BFS::PostIntegration(Snapshot* snapshot, const CVManager& cvmanager)
	{
        auto cvs = cvmanager.GetCVs(cvmask_);        
        std::vector<double> x(cvs.size(),0);

        /*The binned cv space is updated at every step
         *After a certain number of steps has been passed, the system updates a
         *bias projection based on the visited histogram states
         */
        for(size_t i = 0; i < cvs.size(); ++i)
        {
            x[i] = cvs[i]->GetValue();
        }
       
        // The histogram is updated based on the index
        h_->at(x) += 1;
    
        // Update the basis projection after a predefined number of steps
        if(snapshot->GetIteration()  % cyclefreq_ == 0) {
            double beta;
            beta = 1.0 / (snapshot->GetTemperature() * snapshot->GetKb());

            // For systems with poorly defined temperature (ie: 1 particle) the
            // user needs to define their own temperature. This is a hack that
            // will be removed in future versions.

            if(!snapshot->GetTemperature())
            {
                beta = temperature_;
                if(temperature_ == 0)
                {
                    std::cout<<std::endl;
                    std::cerr<<"ERROR: Input temperature needs to be defined for this simulation"<<std::endl;
                    MPI_Abort(world_, EXIT_FAILURE);
                }
            }
            iteration_+= 1;
            ProjectBias(cvs,beta);
            std::cout<<"Node: ["<<mpiid_<<"]"<<std::setw(10)<<"\tSweep: "<<iteration_<<std::endl;
        }

		// This gets the bias force from the grid
        std::vector<double> bias_grad(cvs.size(),0);

        InBounds(cvs);
        if (bounds_) 
        {
            auto& bias = f_->at(x);	
            for (size_t ii = 0; ii<bias.size(); ii++)
                bias_grad[ii] = bias[ii];
        }
        else 
        {
            // This is where the wall potentials are going to be thrown into the method if the system is not a periodic CV
            for(size_t j = 0; j < cvs.size(); ++j)
            {
                if(!h_->GetPeriodic(j))
                {
                    if(x[j] > boundUp_[j])
                        bias_grad[j] = restraint_[j] * (x[j] - boundUp_[j]);
                    else if(x[j] < boundLow_[j])
                        bias_grad[j] = -restraint_[j] * (x[j] - boundLow_[j]);
                }
            }
        }

		// Take each CV and add its biased forces to the atoms using the chain rule
		auto& forces = snapshot->GetForces();
        auto& virial = snapshot->GetVirial();
		for(size_t i = 0; i < cvs.size(); ++i)
		{
			auto& grad = cvs[i]->GetGradient();
            auto& boxgrad = cvs[i]->GetBoxGradient();

			/* Update the forces in snapshot by adding in the force bias from each
			 *CV to each atom based on the gradient of the CV.
             */
			for (size_t j = 0; j < forces.size(); ++j) 
				forces[j] += bias_grad[i]*grad[j];
            
            virial -= bias_grad[i]*boxgrad;
		}
	}

	// Post-simulation hook.
	void BFS::PostSimulation(Snapshot*, const CVManager&)
	{
	    std::cout<<"Run has finished"<<std::endl;	
	}
 
	// Update the coefficients/bias projection
	void BFS::ProjectBias(const CVList& cvs, const double beta)
	{
        double sum  = 0.0;

        // For multiple walkers, the struct is unpacked
        Grid<uint> histlocal(*h_);

        // Summed between all walkers
        MPI_Allreduce(histlocal.data(), h_->data(), h_->size(), MPI_INT, MPI_SUM, world_);

        // Construct the biased histogram
        size_t i = 0;
        for (Grid<uint>::iterator it2 = h_->begin(); it2 != h_->end(); ++it2, ++i)
        {
            // This is to make sure that the CV projects across the entire surface
            if (*it2 == 0) { *it2 = 1; }
            
            /* The evaluation of the biased histogram which projects the histogram to the
             * current bias of CV space.
             */
            unbias_[i] += (*it2) * exp(b_->at(it2.coordinates())) * weight_ / (double)(cyclefreq_);
        }

        // Reset histogram
        for (auto &val : *h_) {
            val = 0;
        }

        // Create the log array and send that to the integrator
        std::vector<double> z_array (unbias_.size(),0);
        for (i = 0; i < unbias_.size(); i++)
            z_array[i] = log(unbias_[i]);

        //Update the coefficients and determine the difference from the previous iteration
        sum = evaluator_.UpdateCoeff(z_array,*h_);
        coeff_arr_ = evaluator_.GetCoeff();
        //Update both the gradient and the bias on the grids
        evaluator_.UpdateBias(b_,f_);

        if(world_.rank() == 0)
            // Write coeff at this step, but only one walker
            PrintBias(cvs,beta);

        // The convergence tolerance and whether the user wants to exit are incorporated here
        if(sum < tol_)
        {
            std::cout<<"System has converged"<<std::endl;
            if(converge_exit_)
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
    void BFS::PrintBias(const CVList& cvs, const double beta)
    {
        // The filenames will have a standard name, with a user-defined suffix
        std::string filename1 = "basis"+bnme_+".out"; 
        std::string filename2 = "coeff"+bnme_+".out"; 
        std::ofstream basisout;
        std::ofstream coeffout;
        basisout.precision(5);
        coeffout.precision(5);
        basisout.open(filename1.c_str());
        coeffout.open(filename2.c_str());

        // The CV values, PMF projection, PMF, and biased histogram are output for the user
        basisout << "CV Values" << std::setw(35*cvs.size()) << "Basis Set Bias" << std::setw(35) << "PMF Estimate" << std::setw(35) << "Biased Histogram" << std::endl;
        coeffout << "Index" << std::setw(35)<< "Coefficients" << std::endl;
        
        size_t j = 0;
        for(Grid<double>::iterator it = b_->begin(); it != b_->end(); ++it, ++j)
        {
            for(size_t k = 0; k < cvs.size(); ++k)
            {
                // Evaluate the CV values for printing purposes
                basisout << it.coordinate(k) << std::setw(35);
            }
            basisout << -(*it) << std::setw(35);
            basisout << -log(unbias_[j]) << std::setw(35);
            basisout << unbias_[j];
            basisout << std::endl;
        }

        std::vector<double> coeff = evaluator_.GetCoeff();
        size_t ii = 0;
        for (auto& val : coeff) {
            coeffout << ii << std::setw(35) << val << std::endl;
            ii++;
        }

		basisout << std::endl;
		coeffout << std::endl;
        basisout.close();
        coeffout.close();
	}

    // The forces are calculated by chain rule, first  the derivatives of the basis set, then in the PostIntegration function, the derivative of the CV is evaluated
	void BFS::InBounds(const CVList& cvs)
	{	
        std::vector<double> x (cvs.size(),0);
        //This is calculating the derivatives for the bias force
        for (size_t j = 0; j < cvs.size(); ++j)
        {
            x[j] = cvs[j]->GetValue();
            double min = h_->GetLower(j);
            double max = h_->GetUpper(j);

            if(!h_->GetPeriodic(j))
            {
                // In order to prevent the index for the histogram from going out of bounds a check is in place
                if(x[j] > max && bounds_)
                {
                    std::cout<<"WARNING: CV is above the maximum boundary."<<std::endl;
                    std::cout<<"Statistics will not be gathered during this interval"<<std::endl;
                    bounds_ = false;
                }
                else if(x[j] < min && bounds_)
                {
                    std::cout<<"WARNING: CV is below the minimum boundary."<<std::endl;
                    std::cout<<"Statistics will not be gathered during this interval"<<std::endl;
                    bounds_ = false;
                }
                else if(x[j] < max && x[j] > min && !bounds_)
                {
                    std::cout<<"CV has returned in between bounds. Run is resuming"<<std::endl;
                    bounds_ = true;
                }
            }
        }
    }

	//! \copydoc Method::Build()
	BFS* BFS::Build(const Json::Value& json, 
			       		    const MPI_Comm& world,
					        const MPI_Comm& comm,
					        const std::string& path)
    {
		ObjectRequirement validator;
		Value schema;
		Reader reader;
        
		reader.parse(JsonSchema::BFSMethod, schema);
        validator.Parse(schema, path);

        //Validate Inputs
        validator.Validate(json, path);
        if(validator.HasErrors())
            throw BuildException(validator.GetErrors());

        std::vector<double> restrCV(0);
        for(auto& restr : json["CV_restraint_spring_constants"])
            restrCV.push_back(restr.asDouble());

        std::vector<double> boundLow(0);
        for(auto& bndl : json["CV_restraint_minimums"])
            boundLow.push_back(bndl.asDouble());
    
        std::vector<double> boundUp(0);
        for(auto& bndu : json["CV_restraint_maximums"])
            boundUp.push_back(bndu.asDouble());

        auto cyclefreq = json.get("cycle_frequency", 100000).asInt();
        auto freq = json.get("frequency", 1).asInt();
        auto wght = json.get("weight", 1.0).asDouble();
        auto bnme = json.get("basis_filename", "").asString();
        auto temp = json.get("temperature", 0.0).asDouble();
        auto tol  = json.get("tolerance", 1e-6).asDouble();
        auto conv = json.get("convergence_exit", false).asBool();

        Grid<uint> *h = Grid<uint>::BuildGrid(
                                        json.get("grid", Json::Value()) );

        Grid<std::vector<double>> *f = Grid<std::vector<double>>::BuildGrid(
                                        json.get("grid", Json::Value()) );

        Grid<double> *b = Grid<double>::BuildGrid(
                                        json.get("grid", Json::Value()) );

        size_t ii = 0;
        std::vector<BasisFunction*> functions;
        for(auto& m : json["basis_functions"])
        {
            auto *bf = BasisFunction::Build(m, path, b->GetNumPoints(ii));
            functions.push_back(bf);
            ii++;
        }

        auto* m = new BFS(world, comm, h, f, b, functions, restrCV, boundUp, boundLow,
                            cyclefreq, freq, bnme, temp, tol, wght,
                            conv);
      
        if(json.isMember("iteration"))
            m->SetIteration(json.get("iteration",0).asInt());

        if(json.isMember("coefficients") && json.isMember("bias_hist"))
        {
            std::vector<double> coeff;
            std::vector<double> unbias;

            for(auto& c : json["coefficients"])
                coeff.push_back(c.asDouble());

            for(auto& u : json["bias_hist"])
                unbias.push_back(u.asDouble());

            m->SetBasis(coeff, unbias);
    	}

		return m;
    }

}

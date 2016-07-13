#pragma once

#include "Method.h"
#include "../CVs/CollectiveVariable.h"
#include "../Grids/Grid.h"
#include <fstream>
#include <vector>


namespace SSAGES
{

    //! Basis Function Sampling Algorithm
    /*!
     * \ingroup Methods
     *
     *
     * Implementation of the Basis Function Sampling Method based on:
     * Whitmer, J.; Chiu, C.; Joshi, A.; de Pablo J. J. Basis Function Sampling: A New Paradigm for Material Property Computation, Phys. Rev. Letters. 2014, 113, 190602
     */
    /*! A clean mapping structure for both the histogram and the coefficients.
     * All vectors are written as 1D with a row major mapping. In order to make iterating easier, the mapping of the 1D vectors are written here.
     */
    struct Map
    {
        //The coefficient value
        double value;

        //The mapping in an array of integers
        std::vector<int> map;

        Map(const std::vector<int>& map,
            double value) :
            map(map), value(value)
        {}
    };

	/*! The structure that holds the Look-up table for the basis function.
     * To prevent repeated calculations, both the derivatives and values of the Legendre polynomials is stored here. More will be added in future versions
     */
	struct BasisLUT
	{
		//The values of the basis sets		
		std::vector<double> values;

		//The values of the derivatives of the basis sets
		std::vector<double> derivs;

		BasisLUT(const std::vector<double>& values,
			const std::vector<double>& derivs) :
			values(values), derivs(derivs)
		{}
	};
		
	class Basis : public Method
	{
	private:	
        
        //! Histogram of visited states that is stored locally. It is a 1D vector that holds N dimensional data over the number of walkers using a row major mapping
        std::vector<Map> _hist;

        //! Locally defined histogram array for easy mpi operations (it does take up more memory)
        std::vector<int> _histlocal;

        //! Globally defined histogram array for easy mpi operations (lots o memory)
        std::vector<int> _histglobal;

		/*! The globally located coefficient values. As coefficients are updated at the same time, the coefficients should only be defined globally
         */
		std::vector<Map> _coeff;

        /*! The biased histogram of states of the form _hist*exp(phi*beta) where phi is the bias potential and beta is the inverse of the temperature.
         * It is defined globally.
         */
        std::vector<double> _unbias;

        //! The Basis set lookup table, also defined globally
		std::vector<BasisLUT> _LUT;

        /*! The derivatives of the bias potential imposed on the system.
         * They are evaluated by using the lookup table
         */
		std::vector<double> _derivatives;

        /*! The order of the basis polynomials
         */
        std::vector<unsigned int> _polyords;

        //! Storing number of bins for simplicity and readability of code
        std::vector<unsigned int> _nbins;

        /*! The system uses this to determine if the system is to be restrained on the defined interval.
         * The user inputs the spring constants if the system is not periodic.
         */
        std::vector<double> _restraint;

        /*! The positions of the string restraints are placed here in both the upper and lower vectors
         */
        std::vector<double> _boundUp;
        std::vector<double> _boundLow;
        
        //! Frequency of coefficient updates
		unsigned int _cyclefreq;

        //! This is a parameter for allowing the read coefficient function to work correctly
        unsigned int _iter;
        
        //! The node that the current system belongs to, primarily for printing and debugging
        unsigned int _mpiid;

        //! Weighting will be used to potentially sample faster, however it can cause the simulation to explode. By default this value will be set to 1.0
        double _weight;

        //! In the case of the MD engine using a poorly defined temperature, the option to throw it into the method is available. If not provided it takes the value from the engine.
        double _temperature;

        //! The tolerance criteria for the system 
        double _tol;

        //! The checker to see if a run is continuing from a previous one 
        bool _read;

        //! A variable to check to see if the simulation is in bounds or not
        bool _bounds;

        //! A check to see if you want the system to end when it reaches the convergence criteria
        bool _converge_exit;

		//! Updates the bias projection of the PMF
		void UpdateBias(const CVList& cvs, double);

		//! Computes the bias force.
		void CalcBiasForce(const CVList& cvs);

		//! Prints the current bias to a defined file from the JSON
		void PrintBias(const CVList& cvs);

        //! Initializes the look up tables for polynomials
        void BasisInit(const CVList& cvs);

        //! Reads the bias files to continue a previous simulation
        void ReadBasis(std::string, std::string);

        //! Determines the quadrature rule for the basis set. Currently only for Legendre polynomials
        void GaussQuad(std::vector<double>&, std::vector<double>&, int, int);

        //! Returns the bins from the grid
        std::vector<int> GetBin(std::vector<float>&);
        
		//! Output stream for basis projection data.
		std::ofstream _basisout;

        //! Output stream for coefficients (for reading purposes)
        std::ofstream _coeffout;

        //! The option to name both the basis and coefficient files will be given
        //! Basis filename 
        std::string _bnme;

        //! Coefficient filename
        std::string _cnme;


	public: 
		/*! Constructs an instance of the Basis function sampling  method.
         *The coeff describes the basis projection of the system. This is
         *updated once every _cyclefreq. For now, only the legendre polynomial
         *is implemented. Others will be added later
         */
		Basis(boost::mpi::communicator& world,
			 boost::mpi::communicator& comm,
			 const std::vector<unsigned int>& polyord,
             const std::vector<double>& restraint,
             const std::vector<double>& boundUp,
             const std::vector<double>& boundLow,
             unsigned int cyclefreq,
			 unsigned int frequency,
             const std::string bnme,
             const std::string cnme,
             const double temperature,
             const double tol,
             const double weight,
             bool read,
             bool converge) : 
		Method(frequency, world, comm), _coeff(), _hist(), _unbias(), _mpiid(0),
		 _derivatives(), _histlocal(), _cyclefreq(cyclefreq), _polyords(polyord),
         _read(read), _weight(weight), _LUT(), _restraint(restraint), _bnme(bnme),
         _cnme(cnme), _temperature(temperature), _histglobal(), _nbins(), _boundUp(boundUp),
         _boundLow(boundLow), _tol(tol), _converge_exit(converge) 
		{
		}

		// Pre-simulation hook.
		void PreSimulation(Snapshot* snapshot, const CVList& cvs) override;

		// Post-integration hook.
		void PostIntegration(Snapshot* snapshot, const CVList& cvs) override;

		// Post-simulation hook.
		void PostSimulation(Snapshot* snapshot, const CVList& cvs) override;

		void Serialize(Json::Value& json) const override
		{

		}

		~Basis() {}
	};
}
			

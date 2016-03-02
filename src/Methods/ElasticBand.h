#pragma once

#include "Method.h"
#include "../CVs/CollectiveVariable.h"
#include <fstream>

namespace SSAGES
{
	// Implementation of a multi-walker Elastic Band
	// method with no bells and whistles.
	class ElasticBand : public Method
	{
	private:	
		
		// Number Iterations steps, number of iterations 
		// of the elastic band method
		unsigned int _iterations;

		// Number Equilibration steps, number of MD steps to
		// allow the system to reequilibrate before evolving
		unsigned int _equilibrate;

		// Number evolution steps, number of MD steps before
		// collecting statistics for gradients 
		unsigned int _evolution;

		// Number samples, number of samples to 
		// average statistics for gradients 
		unsigned int _nsamples;

		//Number samples actually sampled
		unsigned int _nsampled;

		// Gradiant values
		std::vector<double> _gradient;

		// CV starting location values
		std::vector<double> _centers;

		// Vector of spring constants.
		std::vector<double> _kspring;

		// String spring
		double _stringspring;

		// The node this belongs to
		int _mpiid;

		//current field, or center value
		std::vector<double> _curr_field;

		//Time step of string change
		double _timestep;

		// Output stream for string data.
		std::ofstream _stringout;

		// Current iteration of elastic band method
		unsigned int _currentiter;

		// Restart frequency
		unsigned int _restartiter;
		// Adds a new hill.
		void StringUpdate();

		// Prints the new hill to file
		void PrintString(const CVList& CV);

	public: 
		// Constructs an instance of Elastic Band method.
		// Force frequency to be 1? QUESTION
		// Generation of initial string? QUESTION
		// Always umbrella in string method? QUESTION
		// isteps = Number Iterations steps
		// eqsteps = Number Equilibration steps
		// evsteps = Number evolution steps
		ElasticBand(boost::mpi::communicator& world,
					boost::mpi::communicator& comm,
					unsigned int isteps,
					unsigned int eqsteps,
					unsigned int evsteps,
					unsigned int nsamples,
					const std::vector<double>& centers,
					const std::vector<double>& kspring,
					double stringspring,
					double timestep,
			 		unsigned int frequency) : 
		Method(frequency, world, comm), _iterations(isteps), _equilibrate(eqsteps),
		_evolution(evsteps), _nsamples(nsamples), _nsampled(0),
		_gradient(), _centers(centers), _kspring(kspring), 
		_stringspring(stringspring), _mpiid(0), _curr_field(), 
		_timestep(timestep)
		{
		}

		// Pre-simulation hook.
		void PreSimulation(Snapshot* snapshot, const CVList& cvs) override;

		// Post-integration hook.
		void PostIntegration(Snapshot* snapshot, const CVList& cvs) override;

		// Post-simulation hook.
		void PostSimulation(Snapshot* snapshot, const CVList& cvs) override;

		~ElasticBand() {}
	};
}
			

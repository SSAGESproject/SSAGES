#pragma once

#include "Method.h"
#include "../CVs/CollectiveVariable.h"
#include <fstream>

namespace SSAGES
{
	//! Multi-walker Elastic Band.
	/*!
	 * Implementation of a multi-walker Elastic Band method with no bells and
	 * whistles.
	 *
	 * \ingroup Methods
	 */
	class ElasticBand : public Method
	{
	private:	
		
		//! Number Iterations steps, number of iterations
		//! of the elastic band method.
		unsigned int _iterations;

		//! Number Equilibration steps, number of MD steps to
		//! allow the system to reequilibrate before evolving.
		unsigned int _equilibrate;

		//! Number evolution steps, number of MD steps before
		//! collecting statistics for gradients.
		unsigned int _evolution;

		//! Number samples, number of samples to
		//! average statistics for gradients.
		unsigned int _nsamples;

		//! Number samples actually sampled.
		unsigned int _nsampled;

		//! Gradiant values.
		std::vector<double> _gradient;

		//! CV starting location values.
		std::vector<double> _centers;

		//! Vector of spring constants.
		std::vector<double> _kspring;

		//! String spring constant.
		double _stringspring;

		//! The node this belongs to
		int _mpiid;

		//! current field, or center value
		std::vector<double> _curr_field;

		//! Time step of string change
		double _timestep;

		//! Output stream for string data.
		std::ofstream _stringout;

		//! Current iteration of elastic band method
		unsigned int _currentiter;

		//! Restart frequency
		unsigned int _restartiter;

		//! Adds a new hill.
		void StringUpdate();

		//! Prints the new hill to file
		/*!
		 * \param CV List of CVs.
		 */
		void PrintString(const CVList& CV);

	public:
		//! Constructor
		/*!
		 * \param world MPI global communicator.
		 * \param comm MPI local communicator.
		 * \param isteps Number of iteration steps.
		 * \param eqsteps Number of equilibration steps.
		 * \param evsteps Number of evolution steps.
		 * \param nsamples Number of samples.
		 * \param centers Values defining the center.
		 * \param kspring List of spring constants.
		 * \param stringspring Spring constant of elastic band.
		 * \param timestep Simulation timestep.
		 * \param frequency Frequency with which this Method is invoked.
		 *
		 * Constructs an instance of Elastic Band method.
		 *
		 * Force frequency to be 1? QUESTION
		 * Generation of initial string? QUESTION
		 * Always umbrella in string method? QUESTION
		 */
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

		//! Pre-simulation hook.
		/*!
		 * \param snapshot Simulation snapshot.
		 * \param cvs List of CVs.
		 */
		void PreSimulation(Snapshot* snapshot, const CVList& cvs) override;

		//! Post-integration hook.
		/*!
		 * \param snapshot Current simulation snapshot.
		 * \param cvs List of CVs.
		 */
		void PostIntegration(Snapshot* snapshot, const CVList& cvs) override;

		//! Post-simulation hook.
		/*!
		 * \param snapshot Current simulation snapshot.
		 * \param cvs List of CVs.
		 */
		void PostSimulation(Snapshot* snapshot, const CVList& cvs) override;

		//! \copydoc Serializable::Serialize()
		/*!
		 * \warning Serialization not implemented yet!
		 */
		void Serialize(Json::Value& json) const override
		{

		}

		//! Destructor.
		~ElasticBand() {}
	};
}
			

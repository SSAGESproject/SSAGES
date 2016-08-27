#pragma once

#include "StringMethod.h"
#include "../CVs/CollectiveVariable.h"
#include <fstream>

namespace SSAGES
{
	//! Finite Temperature Spring Method
	/*!
	 * \ingroup Methods
	 *
	 * Implementation of a multi-walker finite string
	 * method with hard wall voronoi cells and running block averages.
	 */
	class FiniteTempString : public StringMethod
	{
	private:	

		//! String modification parameter
		double _kappa;

		//! Number of steps to block average the CV's postions over
		unsigned int _blockiterations;

		//! Time step of string change
		double _tau;

		//! Minimum number of steps to apply umbrella sampling.
		unsigned int _min_num_umbrella_steps;

		//! Flag to run umbrella or not during post-integration        
        bool _run_umbrella;

        //! Iterator that keeps track of umbrella iterations
		unsigned int _umbrella_iter;

		//! Stores the last positions of the CVs
		std::vector<double> _prev_CVs;

		//! Steores the last positions of the atoms
		std::vector<Vector3> _prev_positions;

		//! Checks if CV is in voronoi cell
		bool InCell(const CVList& cvs) const;

		//! Updates the string according to the FTS method
		void StringUpdate();

	public:
		//! Constructor
		/*!
		 * \param world MPI global communicator.
		 * \param comm MPI local communicator.
		 * \param centers List of centers.
		 * \param maxiterator Maximum number of iterations.
		 * \param blockiterations Number of iterations per block averaging.
		 * \param tau Value of tau (default: 0.1).
		 * \param cvspring Spring constants for cvs.
		 * \param kappa Value of kappa (default: 0.1).
		 * \param frequency Frequency with which this method is invoked.
		 *
		 * Constructs an instance of Finite String method.
		 */
		FiniteTempString(boost::mpi::communicator& world,
					boost::mpi::communicator& comm,
					const std::vector<double>& centers,
					unsigned int maxiterations,
					unsigned int blockiterations,
					double tau,
					const std::vector<double> cvspring,
					double kappa,
					unsigned int springiter,
			 		unsigned int frequency) : 
		StringMethod(world, comm, centers, maxiterations, cvspring, frequency),
		_kappa(kappa), _blockiterations(blockiterations), _tau(tau), 
		_min_num_umbrella_steps(springiter), _run_umbrella(true),
		_umbrella_iter(1)
        {
		}

		//! Post-integration hook.
		void PostIntegration(Snapshot* snapshot, const CVList& cvs) override;
        
		void Serialize(Json::Value& json) const override
        {
        	StringMethod::Serialize(json);

            json["umbrella_iterations"] = _min_num_umbrella_steps;
            json["flavor"] = "FTS";
            json["kappa"] = _kappa;
            json["block_iterations"] = _blockiterations;
            json["time_step"] = _tau;

            for(auto& nw : _newcenters)
            	json["running_average"].append(nw);
        }

		//! Destructor
		~FiniteTempString() {}
	};
}
			

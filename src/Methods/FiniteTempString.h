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

		unsigned int _spring_iter;
        
        bool _run_SMD;

		bool InCell(const CVList& cvs) const;

		void StringUpdate();

	public:
		//! Constructor
		/*!
		 * \param world MPI global communicator.
		 * \param comm MPI local communicator.
		 * \param centers List of centers.
		 * \param NumNodes Number of nodes.
		 * \param maxiterator Maximum number of iterations.
		 * \param blockiterations Number of iterations per block averaging.
		 * \param tau Value of tau (default: 0.1).
		 * \param cvspring Spring constants for cvs.
		 * \param run_SMD Run steered MD to direct CV to proper starting configuration.
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
			 		unsigned int frequency) : 
		StringMethod(world, comm, centers, maxiterations, cvspring, frequency), _blockiterations(blockiterations), _tau(tau), _kappa(kappa), _spring_iter(1), _run_SMD(false)
        {
			//Force _run_SMD for now. Future release will include more details. 
			_run_SMD = true;
		}

		//! Post-integration hook.
		void PostIntegration(Snapshot* snapshot, const CVList& cvs) override;
        
		void Serialize(Json::Value& json) const override
        {
        	StringMethod::Serialize(json);

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
			

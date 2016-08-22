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

		unsigned int _iterator;

		unsigned int _spring_iter;

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
		StringMethod(world, comm, centers, maxiterations, blockiterations,
		tau, cvspring, frequency), _kappa(kappa), _iterator(1), _spring_iter(1)
		{
		}

		//! Post-integration hook.
		void PostIntegration(Snapshot* snapshot, const CVList& cvs) override;
        
		void Serialize(Json::Value& json) const override
        {
        	StringMethod::Serialize(json);

            json["type"] = "FTS";
            json["kappa"] = _kappa;

            // For now must alwasy be true.
            // This ensures you dont have to print _prev_positions, and _cv_prev
            json["run_SMD"] = true;

        }

		//! Destructor
		~FiniteTempString() {}
	};
}
			

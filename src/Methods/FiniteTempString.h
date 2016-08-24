/**
 * This file is part of
 * SSAGES - Suite for Advanced Generalized Ensemble Simulations
 *
 * Copyright 2016 Ben Sikora <bsikora906@gmail.com>
 *                Cody Bezik <bezik@uchicago.edu>
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
		tau, cvspring, frequency), _kappa(kappa), _spring_iter(1)
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
			

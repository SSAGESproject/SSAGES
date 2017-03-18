/**
 * This file is part of
 * SSAGES - Suite for Advanced Generalized Ensemble Simulations
 *
 * Copyright 2016 Hythem Sidky <hsidky@nd.edu>
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

#include "Driver.h"
#include "../../hooks/qbox/QboxHook.h"

namespace SSAGES
{
	//! Driver for QBox simulation. 
	class QboxDriver : public Driver
	{
	private:
		//! Number of MD steps.
		int mdsteps_;

		//! Number of QM (electronic) iterations. 
		int qmsteps_; 
		
		//! Pointer to QBox hook.
		QboxHook* qbhook_;

	public:

		//! Constructor. 
		/*! 
		 * \param world_comm MPI global communicator.
		 * \param local_comm MPI local communicator.
		 * \param walkerID ID of the walker assigned to this driver.
		 */
		QboxDriver(mpi::communicator& world_comm, 
		           mpi::communicator& local_comm,
		           int walkerid) : 
		Driver(world_comm, local_comm, walkerid)
		{}
	
		//! Run simulation.
		void Run() override;

		//! Run input file.
		void ExecuteInputFile(std::string) override
		{
		}

		//! Set up the driver
		/*!
		 * \param json JSON value containing input information.
		 * \param path Path for JSON path specification.
		 */
		void BuildDriver(const Json::Value& json, const std::string& path) override
		{
			inputfile_ = json.get("inputfile","none").asString();
			mdsteps_ = json.get("MDSteps", 1).asInt();

			// Set hook. 
			qbhook_ = new QboxHook();
			qmsteps_ = json.get("qm_steps", 30).asDouble();
			hook_ = dynamic_cast<Hook*>(qbhook_);
		}

	};
}
/*
 * This file is part of
 * SSAGES - Suite for Advanced Generalized Ensemble Simulations
 *
 * Copyright 2017 Hythem Sidky <hsidky@nd.edu>
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

#include "config.h"
#include EXTRA_CONFIG

#include <boost/any.hpp>

#include "utils/simError.h"
#include "brains/Register.hpp"
#include "brains/SimCreator.hpp"
#include "brains/SimInfo.hpp"
#include "integrators/IntegratorFactory.hpp"
#include "integrators/Integrator.hpp"

#include "Driver.h"
#include "../../hooks/openmd/OpenMDHook.h"

namespace SSAGES
{
	//! Driver for OpenMD simulations.
	class OpenMDDriver : public Driver
	{
	private:
		//! Number of walkers. 
		int nwalkers_; 

		//! Number of MD steps. 
		int MDSteps_; 

		//! Is the simulation a restart? 
		bool restart_;

	public:
		//! Constructor. 
		/*!
		 * \param world_comm MPI global communicator. 
		 * \param local_comm MPI local communicator. 
		 * \param walkerID ID of the walker assigned to this driver. 
		 */
		OpenMDDriver(mpi::communicator& world_comm,
		             mpi::communicator& local_comm,
		             int walkerID) : 
		Driver(world_comm, local_comm, walkerID), 
		nwalkers_(world_comm.size()/local_comm.size())
		{}

		//! Run simulation
		void Run() override
		{
	  		errorCheckPoint();
	  		OpenMD::registerAll();

	  		worldRank = _world.rank();

	  		OpenMD::SimCreator creator;
	  		auto* info = creator.createSim(_inputfile);
	  		auto* simparams = info->getSimParams();
	  		if(!simparams->haveEnsemble())
	  			throw std::invalid_argument("SSAGES only supports ensemble mode for OpenMD.");

 		   auto* integrator = OpenMD::IntegratorFactory::getInstance()->createIntegrator(OpenMD::toUpperCopy(simparams->getEnsemble()), info);
 		   integrator->integrate();
 		   delete integrator;
		}

		//! Run Input file.
		/*!
		 * In the case of OpenMD, the input file does not need to be parsed for
		 * special information.
		 */
		void ExecuteInputFile(std::string) override
		{
		}

		//! Set up the driver.
		/*!
		 * \param json JSON value containing input information.
		 * \param path Path for JSON path specification.
		 */
		void BuildDriver(const Json::Value& json, const std::string& path) override
		{
			_inputfile = json.get("inputfile","none").asString();
			MDSteps_ = json.get("MDSteps", 1).asInt();

			// Set hook. 
			auto& hook = OpenMDHook::Instance();
			//hook.SetIterationTarget(MDSteps_);
			_hook = dynamic_cast<Hook*>(&hook);

			// Restart?
			restart_ = json.get("restart", false).asBool();
		}

		//! \copydoc Serializable::Serialize()
		void Serialize(Json::Value& json) const override
		{
			// Call parent first.
			Driver::Serialize(json);

			json["MDSteps"] = MDSteps_;
			json["type"] = "OpenMD";
		}

	};
}
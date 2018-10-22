/**
 * This file is part of
 * SSAGES - Software Suite for Advanced General Ensemble Simulations
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

#include "Driver.h"
#include "ResourceHandler.h"
#include "OpenMDHook.h"
#include "json/json.h"

#include "config.h"
#include EXTRA_CONFIG
#include "utils/simError.h"
#include "brains/Register.hpp"
#include "brains/SimCreator.hpp"
#include "brains/SimInfo.hpp"
#include "integrators/IntegratorFactory.hpp"
#include "integrators/Integrator.hpp"

namespace SSAGES
{
	void Driver::Run()
	{
		errorCheckPoint();
		OpenMD::registerAll();

		worldRank = rh_->GetWorldMxxComm().rank();

		OpenMD::SimCreator creator;
		auto* info = creator.createSim(rh_->GetInput());
		auto* simparams = info->getSimParams();
		if(!simparams->haveEnsemble())
			throw std::invalid_argument("SSAGES only supports ensemble mode for OpenMD.");

		auto* integrator = OpenMD::IntegratorFactory::getInstance()->createIntegrator(OpenMD::toUpperCopy(simparams->getEnsemble()), info);
		integrator->integrate();
		delete integrator;
	}

	Driver* Driver::Build(const Json::Value& json, const MPI_Comm& world)
	{
		auto* rh = ResourceHandler::Build(json, world);
		auto& hook = OpenMDHook::Instance();
		rh->ConfigureHook(dynamic_cast<Hook*>(&hook));
		return new Driver(rh);
	}   

	Driver::~Driver()
	{
		delete rh_;
	} 
}
/**
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
#include "Method.h"
#include "json/json.h"
#include "schema.h"
#include "Validator/ObjectRequirement.h"
#include "Validator/ArrayRequirement.h"
#include "ABF.h"
#include "Umbrella.h"
#include "BasisFunc.h"
#include "ForwardFlux.h"
#include "Meta.h"
#include "StringMethod.h"
#include <stdexcept>

using namespace Json;

namespace SSAGES
{
	Method* Method::BuildMethod(const Json::Value& json, 
		                        const MPI_Comm& world, 
							    const MPI_Comm& comm, 
							    const std::string& path)
	{
		if(json["type"] == "ABF")
			return ABF::Build(json, world, comm, path);
		else if(json["type"] == "Basis")
			return Basis::Build(json, world, comm, path);
		else if(json["type"] == "ForwardFlux")
			return ForwardFlux::Build(json, world, comm, path);
		else if(json["type"] == "Metadynamics")
			return Meta::Build(json, world, comm, path);
		else if(json["type"] == "Umbrella")
			return Umbrella::Build(json, world, comm, path);
		else if(json["type"] == "String")
			return StringMethod::Build(json, world, comm, path);
		else
			throw std::invalid_argument(path + ": Unknown method type specified.");
	}
}


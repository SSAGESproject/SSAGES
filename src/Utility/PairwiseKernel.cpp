/**
 * This file is part of
 * SSAGES - Software Suite for Advanced General Ensemble Simulations
 *
 * Copyright 2017 Hythem Sidky <hsidky@nd.edu>
 *                Michael Quevillon <mquevill@nd.edu>
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

#include <stdexcept>
#include "Validator/ObjectRequirement.h"
#include "Drivers/DriverException.h"
#include "PairwiseKernel.h"
#include "schema.h"

namespace SSAGES
{
	PairwiseKernel* PairwiseKernel::Build(const Json::Value& json, const std::string& path)
	{
		auto type = json.get("type", "none").asString();
		if(type == "gaussian")
			return GaussianPK::Build(json, path);
		else if(type == "rationalswitch")
			return RationalSwitchPK::Build(json, path);
		else
			throw std::invalid_argument("Invalid pairwise kernel type \"" + type + "\".");
	}

	//! Build GaussianPK from JSON value. 
	/*!
	 * \param json JSON value node. 
	 * 
	 * \return Pointer to new GaussianPK.
	 */
	GaussianPK* GaussianPK::Build(const Json::Value& json, const std::string& path)
	{
		Json::ObjectRequirement validator;
		Json::Value schema;
		Json::CharReaderBuilder rbuilder;
		Json::CharReader* reader = rbuilder.newCharReader();

		reader->parse(JsonSchema::GaussianPK.c_str(),
		              JsonSchema::GaussianPK.c_str() + JsonSchema::GaussianPK.size(),
		              &schema, nullptr);
		validator.Parse(schema, path);

		// Validate inputs.
		validator.Validate(json, path);
		if(validator.HasErrors())
			throw BuildException(validator.GetErrors());
		
		return new GaussianPK(
		            json["mu"].asDouble(), 
		            json["sigma"].asDouble()
		        );
	}
	

	//! Build RationalSwitchPK from JSON value. 
	/*!
	 * \param json JSON value node. 
	 * 
	 * \return Pointer to new RationalSwitchPK.
	 */
	RationalSwitchPK* RationalSwitchPK::Build(const Json::Value& json, const std::string& path)
	{
		Json::ObjectRequirement validator;
		Json::Value schema;
		Json::CharReaderBuilder rbuilder;
		Json::CharReader* reader = rbuilder.newCharReader();

		reader->parse(JsonSchema::RationalSwitchPK.c_str(),
		              JsonSchema::RationalSwitchPK.c_str() + JsonSchema::RationalSwitchPK.size(),
		              &schema, nullptr);
		validator.Parse(schema, path);

		// Validate inputs.
		validator.Validate(json, path);
		if(validator.HasErrors())
			throw BuildException(validator.GetErrors());
		return new RationalSwitchPK(
		            json["d0"].asDouble(), 
		            json["r0"].asDouble(), 
		            json["n"].asInt(), 
		            json["m"].asInt()
		        );
	}
}

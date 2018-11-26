/**
 * This file is part of
 * SSAGES - Software Suite for Advanced General Ensemble Simulations
 *
 * Copyright 2016 Jiyuan Li <jyli@uchicago.edu>
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

#include "Constraint.h"
#include "json/json.h"
#include "schema.h"
#include "../Drivers/DriverException.h"
#include "../Validator/ObjectRequirement.h"
#include "../Validator/ArrayRequirement.h"
#include "COPSS.h"
#include "COPSSImage.h"
using namespace Json;

namespace SSAGES
{
	Constraint* Constraint::BuildConstraint(const Json::Value &json,
						boost::mpi::communicator& comm)
	{
		return BuildConstraint(json, comm, "#/Constraints");
	}

	Constraint* Constraint::BuildConstraint(const Json::Value &json,
						boost::mpi::communicator& comm, 
						const std::string& path)
	{
		ObjectRequirement validator;
		Value schema;
		Reader reader;

		Constraint* constraint = nullptr;

		// Random device for seed generation. 
		// std::random_device rd;
		// auto maxi = std::numeric_limits<int>::max();
		// auto seed = json.get("seed", rd() % maxi).asUInt();

		// Get move type. 
		std::string type = json.get("type", "none").asString();

		if(type == "COPSS")
		{
			reader.parse(JsonSchema::COPSSConstraint, schema);
			validator.Parse(schema, path);

			// Validate inputs.
			validator.Validate(json, path);
			if(validator.HasErrors())
				throw BuildException(validator.GetErrors());

			auto* c = new COPSS(comm, 1);

			constraint = static_cast<Constraint*>(c);
		}
		else if(type == "COPSSImage")
		{
			reader.parse(JsonSchema::COPSSImageConstraint, schema);
			validator.Parse(schema, path);
			
			// Validate inputs.
			validator.Validate(json, path);
			if(validator.HasErrors())
				throw BuildException(validator.GetErrors());
			
			double einner = json.get("einner",1).asDouble();

			int ion_type_start = json.get("ion_type_start",1).asInt();
		
			std::vector<double> atomTypeRadius;
			for(auto& atomType : json["atom type radius"])
				atomTypeRadius.push_back(atomType.asDouble());
	
			auto* c = new COPSSImage(comm, 1, einner, ion_type_start, atomTypeRadius);

			constraint = static_cast<Constraint*>(c);

		}
		else
		{
			throw BuildException({path + ": Unknown constraint type specified."+type+" is not a valid type!"});
		}

		return constraint;
	}

	void Constraint::BuildConstraint(const Json::Value &json, 
						  ConstraintList &clist,
						  boost::mpi::communicator& comm,
						  const std::string &path)
	{
		ArrayRequirement validator;
		Value schema;
		Reader reader;

		reader.parse(JsonSchema::constraints, schema);
		validator.Parse(schema, path);

		// Validate high level schema.
		validator.Validate(json, path);
		if(validator.HasErrors())
			throw BuildException(validator.GetErrors());

		// Loop through CVs.
		int i = 0;
		for(auto& m : json)
		{
			clist.push_back(BuildConstraint(m, comm, path + "/" + std::to_string(i)));
			++i;
		}
	}
}

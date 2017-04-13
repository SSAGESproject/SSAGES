/**
 * This file is part of
 * SSAGES - Suite for Advanced Generalized Ensemble Simulations
 *
 * Copyright 2016 Ben Sikora <bsikora906@gmail.com>
 *                Emre Sevgen <sesevgen@uchicago.edu>
 *                Yamil Colon <yamilcolon2015@u.northwestern.edu>
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
#include "CollectiveVariable.h"
#include "AngleCV.h"
#include "BoxVolumeCV.h"
#include "json/json.h"
//#include "ParticleCoordinateCV.h"
//#include "ParticlePositionCV.h"
//#include "ParticleSeparationCV.h"
//#include "TorsionalCV.h"
//#include "GyrationTensorCV.h"
//#include "RMSDCV.h"
//#include "RouseModeCV.h"
//#include "CoordinationNumberCV.h"

namespace SSAGES
{
	CollectiveVariable* CollectiveVariable::BuildCV(const Json::Value &json, const std::string& path)
	{
		// Get move type. 
		auto type = json.get("type", "none").asString();

		if(type == "Angle")
			return AngleCV::Build(json, path);

		/*
		if(type == "ParticleCoordinate")
		{
			reader.parse(JsonSchema::ParticleCoordinateCV, schema);
			validator.Parse(schema, path);

			// Validate inputs.
			validator.Validate(json, path);
			if(validator.HasErrors())
				throw BuildException(validator.GetErrors());

			Label atomids;
			for(auto& id : json["atom_ids"])
				atomids.push_back(id.asInt());

			auto indextype = json.get("dimension","x").asString();

			Dimension index;
			if(indextype == "x")
				index = Dimension::x;
			else if(indextype == "y")
				index = Dimension::y;
			else if(indextype == "z")
				index = Dimension::z;
			else
				throw BuildException({"Could not obtain ParticleCoordinate dimension specified."});

			auto* c = new ParticleCoordinateCV(atomids, index);
			cv = static_cast<CollectiveVariable*>(c);
		}
		else if(type == "ParticlePosition")
		{
			reader.parse(JsonSchema::ParticlePositionCV, schema);
			validator.Parse(schema, path);

			// Validate inputs.
			validator.Validate(json, path);
			if(validator.HasErrors())
				throw BuildException(validator.GetErrors());
			
			Label atomids;
			for(auto& id : json["atom_ids"])
				atomids.push_back(id.asInt());
			
			Vector3 position;
			position[0] = json["position"][0].asDouble();
			position[1] = json["position"][1].asDouble();
			position[2] = json["position"][2].asDouble();

			auto fixx = json["fix"][0].asBool();
			auto fixy = json["fix"][1].asBool();
			auto fixz = json["fix"][2].asBool();

			auto* c = new ParticlePositionCV(atomids, position, fixx, fixy, fixz);

			cv = static_cast<CollectiveVariable*>(c);
		}
		else if(type == "Torsional")
		{
			reader.parse(JsonSchema::TorsionalCV, schema);
			validator.Parse(schema, path);

			// Validate inputs.
			validator.Validate(json, path);
			if(validator.HasErrors())
				throw BuildException(validator.GetErrors());

			std::vector<int> atomids;
			for(auto& s : json["atom_ids"])
				atomids.push_back(s.asInt());

			auto periodic = json.get("periodic", true).asBool();

			auto* c = new TorsionalCV(atomids[0], atomids[1], atomids[2], atomids[3], periodic);

			cv = static_cast<CollectiveVariable*>(c);
		}
		else if(type == "ParticleSeparation")
		{
			reader.parse(JsonSchema::ParticleSeparationCV, schema);
			validator.Parse(schema, path);

			// Validate inputs.
			validator.Validate(json, path);
			if(validator.HasErrors())
				throw BuildException(validator.GetErrors());
			
			std::vector<int> group1, group2;
			
			for(auto& s : json["group1"])
				group1.push_back(s.asInt());

			for(auto& s : json["group2"])
				group2.push_back(s.asInt());

			ParticleSeparationCV* c;
			if(json.isMember("dimension"))
			{
				auto fixx = json["dimension"][0].asBool();
				auto fixy = json["dimension"][1].asBool();
				auto fixz = json["dimension"][2].asBool();

				c = new ParticleSeparationCV(group1, group2, fixx, fixy, fixz);

			}
			else
			{
				c = new ParticleSeparationCV(group1, group2);
			}

			cv = static_cast<CollectiveVariable*>(c);
		}
		else if(type == "GyrationTensor")
		{
			reader.parse(JsonSchema::GyrationTensorCV, schema); 
			validator.Parse(schema, path); 

			// Validate inputs. 
			validator.Validate(json, path); 
			if(validator.HasErrors())
				throw BuildException(validator.GetErrors());

			std::vector<int> atomids; 
			for(auto& s : json["atom_ids"])
				atomids.push_back(s.asInt());

			GyrationTensor component = Rg;
			auto comp = json["component"].asString();
			
			if(comp == "Rg")
				component = Rg; 
			else if(comp == "principal1")
				component = principal1;
			else if(comp == "principal2")
				component = principal2;
			else if(comp == "principal3")
				component = principal3;
			else if(comp == "asphericity")
				component = asphericity;
			else if(comp == "acylindricity")
				component = acylindricity;
			else if(comp == "shapeaniso")
				component = shapeaniso;

			auto* c = new GyrationTensorCV(atomids, component);
			cv = static_cast<CollectiveVariable*>(c);
		}
		else if(type == "RMSD")
		{
			reader.parse(JsonSchema::RMSDCV, schema);
			validator.Parse(schema, path);

			// Validate inputs.
			validator.Validate(json, path);
			if(validator.HasErrors())
				throw BuildException(validator.GetErrors());
			
			std::vector<int> atomids;
			for(auto& s : json["atom ids"])
				atomids.push_back(s.asInt());
			auto reference = json.get("reference"," ").asString(); 

			auto* c = new RMSDCV(atomids, reference, json.get("use_range", false).asBool());
			
			cv = static_cast<CollectiveVariable*>(c);
		}
		else if(type == "RouseMode")
		{
			reader.parse(JsonSchema::RouseModeCV, schema);
			validator.Parse(schema, path);

			// Validate inputs.
			validator.Validate(json, path);
			if(validator.HasErrors())
				throw BuildException(validator.GetErrors());
			
			std::vector<Label> groups;
			for (auto& group : json["groups"]) {
				groups.push_back({});
				for (auto& id : group) {
					groups.back().push_back(id.asInt());
				}
			}
			auto mode = json.get("mode",0).asInt();

			auto* c = new RouseModeCV( groups, mode);
			
			cv = static_cast<CollectiveVariable*>(c);
		}
		else if(type == "CoordinationNumber")
		{
			reader.parse(JsonSchema::CoordinationNumberCV, schema);
			validator.Parse(schema, path);

			// Validate inputs.
			validator.Validate(json, path);
			if(validator.HasErrors())
				throw BuildException(validator.GetErrors());
			
			std::vector<int> group1, group2;
			
			for(auto& s : json["group1"])
				group1.push_back(s.asInt());

			for(auto& s : json["group2"])
				group2.push_back(s.asInt());
			
			auto* c = new CoordinationNumberCV(group1, group2, SwitchingFunction::Build(json["switching"]));
			cv = static_cast<CollectiveVariable*>(c);
		}
		else if(type == "BoxVolume")
		{
			
		}
		else
		{
			throw BuildException({path + ": Unknown CV type specified."+type+" is not a valid type!"});
		}

		*/
	}
}

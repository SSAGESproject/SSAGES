#include "CollectiveVariable.h"
#include "json/json.h"
#include "schema.h"
#include "../Drivers/DriverException.h"
#include "../Validator/ObjectRequirement.h"
#include "../Validator/ArrayRequirement.h"
#include "AtomCoordinateCV.h"
#include "AtomPositionCV.h"
#include "ImproperCV.h"
#include "TorsionalCV.h"

using namespace Json;

namespace SSAGES
{
	CollectiveVariable* CollectiveVariable::BuildCV(const Json::Value &json)
	{
		return BuildCV(json, "#/cvs");
	}

	CollectiveVariable* CollectiveVariable::BuildCV(const Value &json, 
						const std::string& path)
	{
		ObjectRequirement validator;
		Value schema;
		Reader reader;

		CollectiveVariable* cv = nullptr;

		// Random device for seed generation. 
		// std::random_device rd;
		// auto maxi = std::numeric_limits<int>::max();
		// auto seed = json.get("seed", rd() % maxi).asUInt();

		// Get move type. 
		std::string type = json.get("type", "none").asString();

		if(type == "AtomCoordinate")
		{
			reader.parse(JsonSchema::AtomCoordinateCV, schema);
			validator.Parse(schema, path);

			// Validate inputs.
			validator.Validate(json, path);
			if(validator.HasErrors())
				throw BuildException(validator.GetErrors());

			auto atomid = json.get("atom id", -1).asInt();
			auto indextype = json.get("dimension","x").asString();

			int index = -1;

			if(indextype == "x")
				index = 0;
			else if(indextype == "y")
				index = 1;
			else if(indextype == "z")
				index = 2;
			else
			{
				std::cout<<"Error getting dimension for Atom Coordinate CV"<<std::endl;
				exit(0);
			}

			auto* c = new AtomCoordinateCV(atomid, index);

			cv = static_cast<CollectiveVariable*>(c);
		}
		else if(type == "AtomPosition")
		{
			reader.parse(JsonSchema::AtomPositionCV, schema);
			validator.Parse(schema, path);

			// Validate inputs.
			validator.Validate(json, path);
			if(validator.HasErrors())
				throw BuildException(validator.GetErrors());
			
			auto atomid = json.get("atom id", -1).asInt();
			Vector3 position;

			position[0]=json["position"][0].asDouble();
			position[1]=json["position"][1].asDouble();
			position[2]=json["position"][2].asDouble();

			auto fixx = json.get("fixx", false).asBool();
			auto fixy = json.get("fixy", false).asBool();
			auto fixz = json.get("fixz", false).asBool();

			auto* c = new AtomPositionCV(atomid, position, fixx, fixy, fixz);

			cv = static_cast<CollectiveVariable*>(c);
		}
		else if(type == "Improper")
		{
			reader.parse(JsonSchema::ImproperCV, schema);
			validator.Parse(schema, path);

			// Validate inputs.
			validator.Validate(json, path);
			if(validator.HasErrors())
				throw BuildException(validator.GetErrors());

			std::vector<int> atomids;
			for(auto& s : json["atom ids"])
				atomids.push_back(s.asInt());

			auto periodic = json.get("periodic", false).asBool();

			auto* c = new ImproperCV(atomids[0], atomids[1], atomids[2], atomids[3], periodic);

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
			for(auto& s : json["atom ids"])
				atomids.push_back(s.asInt());

			auto periodic = json.get("periodic", false).asBool();

			auto* c = new TorsionalCV(atomids[0], atomids[1], atomids[2], atomids[3], periodic);

			cv = static_cast<CollectiveVariable*>(c);
		}
		else
		{
			throw BuildException({path + ": Unknown method type specified."});
		}

		return cv;
	}
}
#include "Grid.h"
#include "json/json.h"
#include "schema.h"
#include "../Drivers/DriverException.h"
#include "../Validator/ObjectRequirement.h"
#include "../Validator/ArrayRequirement.h"
#include "Grid1D.h"
#include "Grid2D.h"
#include "Grid3D.h"

using namespace Json;

namespace SSAGES
{
	Grid* Grid::BuildGrid(const Json::Value &json)
	{
		return BuildGrid(json, "#/Grids");
	}

	Grid* Grid::BuildGrid(const Value &json, 
						const std::string& path)
	{
		ArrayRequirement validator;
		Value schema;
		Value gridjson;
		Reader reader;

		gridjson = json.get("grid",Json::arrayValue);
		Grid* grid = nullptr;

		if(!json.isMember("grid"))
			return nullptr;

		reader.parse(JsonSchema::grid, schema);
		validator.Parse(schema, path);

		// Validate inputs.
		validator.Validate(gridjson, path);
		if(validator.HasErrors())
			throw BuildException(validator.GetErrors());

		std::vector<double> lower;
		std::vector<double> upper;
		std::vector<bool> periodic;
		std::vector<int> num_points;

		for(auto&m : gridjson)
		{
			lower.push_back(m.get("lower",0.0).asDouble());
			upper.push_back(m.get("upper",0.0).asDouble());
			periodic.push_back(m.get("periodic",false).asBool());
			num_points.push_back(m.get("number points",0).asInt());
		}

		if(lower.size() != upper.size() || lower.size() != periodic.size() || 
			lower.size() != num_points.size())
			throw BuildException({"Grid variables dimensions not the same!"});

		if(lower.size() == 1)
		{
			auto* g = new Grid1D(lower, upper, periodic, num_points);
			grid = static_cast<Grid*>(g);
		}
		else if(lower.size() == 2)
		{
			auto* g = new Grid2D(lower, upper, periodic, num_points);
			grid = static_cast<Grid*>(g);
		}
		else if(lower.size() == 3)
		{
			auto* g = new Grid3D(lower, upper, periodic, num_points);
			grid = static_cast<Grid*>(g);
		}
		else
		{
			//throw build error
		}

		return grid;
	}
}


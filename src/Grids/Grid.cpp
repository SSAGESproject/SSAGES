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

	Grid* Grid::BuildGrid(const Json::Value &json,
						const std::string& path)
	{
		ObjectRequirement validator;
		Value schema;
		Value gridjson;
		Reader reader;

		gridjson = json.get("grid",Json::objectValue);
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

		for(auto&m : gridjson["lower"])
			lower.push_back(m.asDouble());

		for(auto&m : gridjson["upper"])
			upper.push_back(m.asDouble());

		for(auto&m : gridjson["number_points"])
			num_points.push_back(m.asDouble());

		if(lower.size() != upper.size() || lower.size() != num_points.size())
			throw BuildException({"Grid variables dimensions not the same!"});

		if(lower.size() == 1)
		{
			auto* g = new Grid1D(lower, upper, num_points);
			grid = static_cast<Grid*>(g);
		}
		else if(lower.size() == 2)
		{
			auto* g = new Grid2D(lower, upper, num_points);
			grid = static_cast<Grid*>(g);
		}
		else if(lower.size() == 3)
		{
			auto* g = new Grid3D(lower, upper, num_points);
			grid = static_cast<Grid*>(g);
		}
		else
			throw BuildException({"SSAGES currently only accepts 1,2, or 3 dimension grids."});

		if(gridjson.isMember("values"))
		{
			std::vector<double> first_values;
			for(auto& p : gridjson["values"])
				first_values.push_back(p.asDouble());

			grid->SetGrid(first_values);

		}

		if(gridjson.isMember("periodic"))
		{
			for(auto&m : gridjson["periodic"])
				periodic.push_back(m.asBool());

			grid->SetPeriodic(periodic);
		}

		return grid;
	}
}


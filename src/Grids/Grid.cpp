/**
 * This file is part of
 * SSAGES - Suite for Advanced Generalized Ensemble Simulations
 *
 * Copyright 2016 Ben Sikora <bsikora906@gmail.com>
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

#include "Grid.h"
#include "json/json.h"
#include "schema.h"
#include "../Drivers/DriverException.h"
#include "../Validator/ObjectRequirement.h"
#include "../Validator/ArrayRequirement.h"

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

	    throw BuildException({"SSAGES currently doesn't accept grids. Work is in progress."});

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


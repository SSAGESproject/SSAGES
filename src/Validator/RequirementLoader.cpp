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
#include "RequirementLoader.h"
#include "StringRequirement.h"
#include "IntegerRequirement.h"
#include "NumberRequirement.h"
#include "ObjectRequirement.h"
#include "ArrayRequirement.h"
#include "BooleanRequirement.h"
#include "NullRequirement.h"
#include "EnumRequirement.h"
#include "AllOfRequirement.h"
#include "AnyOfRequirement.h"
#include "OneOfRequirement.h"
#include "NotRequirement.h"

namespace Json
{
	std::unique_ptr<Requirement> RequirementLoader::LoadRequirement(const Value &json)
	{
		std::unique_ptr<Requirement> item = nullptr;
		// Load up apitemriate requirement type.
		if(json["type"].asString() == "string")
			item = std::move(std::unique_ptr<Requirement>(new StringRequirement()));
		else if(json["type"].asString() == "integer")
			item = std::move(std::unique_ptr<Requirement>(new IntegerRequirement()));
		else if(json["type"].asString() == "number")
			item = std::move(std::unique_ptr<Requirement>(new NumberRequirement()));
		else if(json["type"].asString() == "object")
			item = std::move(std::unique_ptr<Requirement>(new ObjectRequirement()));
		else if(json["type"].asString() == "array")
			item = std::move(std::unique_ptr<Requirement>(new ArrayRequirement()));
		else if(json["type"].asString() == "boolean")
			item = std::move(std::unique_ptr<Requirement>(new BooleanRequirement()));
		else if(json["type"].asString() == "null")
			item = std::move(std::unique_ptr<Requirement>(new NullRequirement()));
		else if(json["allOf"].isArray())
			item = std::move(std::unique_ptr<Requirement>(new AllOfRequirement()));
		else if(json["anyOf"].isArray())
			item = std::move(std::unique_ptr<Requirement>(new AnyOfRequirement()));
		else if(json["oneOf"].isArray())
			item = std::move(std::unique_ptr<Requirement>(new OneOfRequirement()));
		else if(json["enum"].isArray())
			item = std::move(std::unique_ptr<Requirement>(new EnumRequirement()));
		else if(json["not"].isObject())
			item = std::move(std::unique_ptr<Requirement>(new NotRequirement()));
		// last resort.
		else if(json.isObject())
			item = std::move(std::unique_ptr<Requirement>(new ObjectRequirement()));

		return item;
	}

	std::unique_ptr<Requirement> RequirementLoader::LoadExtended(const Value &json)
	{
		std::unique_ptr<Requirement> item = nullptr;
		if(json.isObject() && json.isMember("allOf") && json["allOf"].isArray())
			item = std::move(std::unique_ptr<Requirement>(new AllOfRequirement()));
		else if(json.isObject() && json.isMember("anyOf") && json["anyOf"].isArray())
			item = std::move(std::unique_ptr<Requirement>(new AnyOfRequirement()));
		else if(json.isObject() && json.isMember("oneOf") && json["oneOf"].isArray())
			item = std::move(std::unique_ptr<Requirement>(new OneOfRequirement()));

		return item;
	}

	std::unique_ptr<Requirement> RequirementLoader::LoadRequirement(const ValueType& type)
	{
		std::unique_ptr<Requirement> item = nullptr; 

		switch(type)
		{
			case nullValue:
				item = std::move(std::unique_ptr<Requirement>(new NullRequirement()));
				break;
			case intValue:
				item = std::move(std::unique_ptr<Requirement>(new IntegerRequirement()));
				break;
			case realValue:
				item = std::move(std::unique_ptr<Requirement>(new NumberRequirement()));
				break;
			case stringValue:
				item = std::move(std::unique_ptr<Requirement>(new StringRequirement()));
				break;
			case booleanValue:
				item = std::move(std::unique_ptr<Requirement>(new BooleanRequirement()));
				break;
			case arrayValue:
				item = std::move(std::unique_ptr<Requirement>(new ArrayRequirement()));
				break;
			case objectValue:
				item = std::move(std::unique_ptr<Requirement>(new ObjectRequirement()));
				break;
			default:
				item = nullptr;
		}

		return item;
	}
}
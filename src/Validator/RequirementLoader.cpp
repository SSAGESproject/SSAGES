/**
 * This file has been obtained from
 * SAPHRON - Statistical Applied PHysics through Random On-the-fly Numerics
 * https://github.com/hsidky/SAPHRON
 *
 * Copyright 2016 Hythem Sidky
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE. 
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
			item = std::unique_ptr<Requirement>(new StringRequirement());
		else if(json["type"].asString() == "integer")
			item = std::unique_ptr<Requirement>(new IntegerRequirement());
		else if(json["type"].asString() == "number")
			item = std::unique_ptr<Requirement>(new NumberRequirement());
		else if(json["type"].asString() == "object")
			item = std::unique_ptr<Requirement>(new ObjectRequirement());
		else if(json["type"].asString() == "array")
			item = std::unique_ptr<Requirement>(new ArrayRequirement());
		else if(json["type"].asString() == "boolean")
			item = std::unique_ptr<Requirement>(new BooleanRequirement());
		else if(json["type"].asString() == "null")
			item = std::unique_ptr<Requirement>(new NullRequirement());
		else if(json["allOf"].isArray())
			item = std::unique_ptr<Requirement>(new AllOfRequirement());
		else if(json["anyOf"].isArray())
			item = std::unique_ptr<Requirement>(new AnyOfRequirement());
		else if(json["oneOf"].isArray())
			item = std::unique_ptr<Requirement>(new OneOfRequirement());
		else if(json["enum"].isArray())
			item = std::unique_ptr<Requirement>(new EnumRequirement());
		else if(json["not"].isObject())
			item = std::unique_ptr<Requirement>(new NotRequirement());
		// last resort.
		else if(json.isObject())
			item = std::unique_ptr<Requirement>(new ObjectRequirement());

		return item;
	}

	std::unique_ptr<Requirement> RequirementLoader::LoadExtended(const Value &json)
	{
		std::unique_ptr<Requirement> item = nullptr;
		if(json.isObject() && json.isMember("allOf") && json["allOf"].isArray())
			item = std::unique_ptr<Requirement>(new AllOfRequirement());
		else if(json.isObject() && json.isMember("anyOf") && json["anyOf"].isArray())
			item = std::unique_ptr<Requirement>(new AnyOfRequirement());
		else if(json.isObject() && json.isMember("oneOf") && json["oneOf"].isArray())
			item = std::unique_ptr<Requirement>(new OneOfRequirement());

		return item;
	}

	std::unique_ptr<Requirement> RequirementLoader::LoadRequirement(const ValueType& type)
	{
		std::unique_ptr<Requirement> item = nullptr; 

		switch(type)
		{
			case nullValue:
				item = std::unique_ptr<Requirement>(new NullRequirement());
				break;
			case intValue:
				item = std::unique_ptr<Requirement>(new IntegerRequirement());
				break;
			case realValue:
				item = std::unique_ptr<Requirement>(new NumberRequirement());
				break;
			case stringValue:
				item = std::unique_ptr<Requirement>(new StringRequirement());
				break;
			case booleanValue:
				item = std::unique_ptr<Requirement>(new BooleanRequirement());
				break;
			case arrayValue:
				item = std::unique_ptr<Requirement>(new ArrayRequirement());
				break;
			case objectValue:
				item = std::unique_ptr<Requirement>(new ObjectRequirement());
				break;
			default:
				item = nullptr;
		}

		return item;
	}
}
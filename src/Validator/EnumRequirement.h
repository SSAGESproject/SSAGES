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
#pragma once 

#include "Requirement.h"

namespace Json
{
	//! Requires entry to be member of an enum.
	/*!
	 * \ingroup Json
	 */
	class EnumRequirement : public Requirement
	{
	private:
		std::vector<Value> _enum; //!< Enum value.
	
	public:
		//! Constructor
		EnumRequirement() : _enum(0) {}

		//! Clear enum value.
		virtual void Reset() 
		{
			_enum.clear();
		}

		//! Parse JSON input value to generate enum.
		/*!
		 * \param json JSON input value.
		 */
		virtual void Parse(Value json, const std::string&) 
		{
			if(json.isArray())
				for(auto& val : json)
					_enum.push_back(val);
		}

		//! Validate that JSON value is member of the parsed enum.
		/*!
		 * \param json JSON value to be validated.
		 * \param path Path for JSON path specification.
		 */
		virtual void Validate(const Value& json, const std::string& path)
		{
			bool found = false;

			for(auto& val : _enum)
				if(json == val)
					found = true;

			if(!found)
				PushError(path  + ": Value is not a valid entry.");
		}
	};
}
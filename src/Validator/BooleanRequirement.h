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
	//! Requires json value to be of type Bool.
	/*!
	 * \ingroup Json
	 */
	class BooleanRequirement : public Requirement
	{
	public:
		//! Reset Requirement.
		virtual void Reset() {}

		//! Parse JSON string
		virtual void Parse(Value, const std::string&) {}

		//! Validate that JSON value is Bool.
		/*!
		 * \param json JSON value to validate.
		 * \param path Path for JSON path specification.
		 */
		virtual void Validate(const Value& json, const std::string& path)
		{
			if(!json.isBool())
				PushError(path + ": Must be a boolean");
		}
	};
}
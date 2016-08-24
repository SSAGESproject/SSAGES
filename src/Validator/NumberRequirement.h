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

#include <cmath>
#include "Requirement.h"

namespace Json
{
	//! Requirements on a numeric value
	/*!
	 * The numbers are stored internally as \c double.
	 *
	 * \ingroup Json
	 */
	class NumberRequirement : public Requirement
	{
	private:
		std::string _path; //!< JSON path.
		double _multipleOf; //!< Base value for "multiple of" requirement.
		double _min; //!< Lower bound for range requirement.
		double _max; //!< Upper bound for range requirement.
		bool _multSet; //!< If \c True, "Multiple of" requirement is active.
		bool _minSet; //!< If \c True, Lower bound for range requirement is active.
		bool _maxSet; //!< If \c True, Upper bound for range requirement is active.
		bool _exclMin; //!< If \c True, lower bound is exclusive.
		bool _exclMax; //!< If \c True, upper bound is exclusive.


	public:
		//! Constructor.
		NumberRequirement() : 
		_path(), _multipleOf(0), _min(0), _max(0), _multSet(false), 
		_minSet(false), _maxSet(false), _exclMin(false), _exclMax(false)
		{}

		//! Reset Requirement.
		virtual void Reset() override
		{
			_multipleOf = 0;
			_minSet = _maxSet = false;
			_exclMin = _exclMax = false; 
			_min = _max = 0;
			_multSet = false;
			ClearErrors();
			ClearNotices();
		}

		//! Parse JSON value to set up Requirement.
		/*!
		 * \param json JSON input value.
		 * \param path Path for JSON path specification.
		 */
		virtual void Parse(Value json, const std::string& path) override
		{
			Reset();
			
			_path = path;
			if(json.isMember("multipleOf") && json["multipleOf"].isNumeric())
			{
				_multSet = true;
				_multipleOf = json["multipleOf"].asDouble();
			}

			if(json.isMember("minimum") && json["minimum"].isNumeric())
			{
				_minSet = true;
				_min = json["minimum"].asDouble();
			}

			if(json.isMember("maximum") && json["maximum"].isNumeric())
			{
				_maxSet = true;
				_max = json["maximum"].asDouble();
			}

			if(json.isMember("exclusiveMinimum") && json["exclusiveMinimum"].isBool())
			{
				_exclMin = json["exclusiveMinimum"].asBool();
			}

			if(json.isMember("exclusiveMaximum") && json["exclusiveMaximum"].isBool())
			{
				_exclMax = json["exclusiveMaximum"].asBool();
			}
		}

		//! Validate JSON value.
		/*!
		 * \param json JSON value to validate.
		 * \param path Path for JSON path specification.
		 *
		 * Test that the JSON value meets the requirements set via
		 * NumberRequirement::Parse(). If the validation fails, an error is added
		 * to the list of error messages.
		 */
		virtual void Validate(const Value& json, const std::string& path) override
		{
			if(!json.isNumeric())
			{
				PushError(path + ": Must be of type \"number\".");
				return;
			}

			if(_multSet && fmod(json.asDouble(), _multipleOf) != 0)
				PushError(path + ": Value must be a multiple of " + std::to_string(_multipleOf));

			if(_minSet)
			{
				if(_exclMin && json.asDouble() <= _min)
					PushError(path + ": Value must be greater than " + std::to_string(_min));
				else if(json.asDouble() < _min)
					PushError(path + ": Value cannot be less than " + std::to_string(_min));
			}

			if(_maxSet)
			{
				if(_exclMax && json.asDouble() >= _max)
					PushError(path + ": Value must be less than " + std::to_string(_max));
				else if(json.asDouble() > _max)
					PushError(path + ": Value cannot be greater than " + std::to_string(_max));
			}
		}
	};
}
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
		std::string path_; //!< JSON path.
		double multipleOf_; //!< Base value for "multiple of" requirement.
		double min_; //!< Lower bound for range requirement.
		double max_; //!< Upper bound for range requirement.
		bool multSet_; //!< If \c True, "Multiple of" requirement is active.
		bool minSet_; //!< If \c True, Lower bound for range requirement is active.
		bool maxSet_; //!< If \c True, Upper bound for range requirement is active.
		bool exclMin_; //!< If \c True, lower bound is exclusive.
		bool exclMax_; //!< If \c True, upper bound is exclusive.


	public:
		//! Constructor.
		NumberRequirement() : 
		path_(), multipleOf_(0), min_(0), max_(0), multSet_(false), 
		minSet_(false), maxSet_(false), exclMin_(false), exclMax_(false)
		{}

		//! Reset Requirement.
		virtual void Reset() override
		{
			multipleOf_ = 0;
			minSet_ = maxSet_ = false;
			exclMin_ = exclMax_ = false; 
			min_ = max_ = 0;
			multSet_ = false;
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
			
			path_ = path;
			if(json.isMember("multipleOf") && json["multipleOf"].isNumeric())
			{
				multSet_ = true;
				multipleOf_ = json["multipleOf"].asDouble();
			}

			if(json.isMember("minimum") && json["minimum"].isNumeric())
			{
				minSet_ = true;
				min_ = json["minimum"].asDouble();
			}

			if(json.isMember("maximum") && json["maximum"].isNumeric())
			{
				maxSet_ = true;
				max_ = json["maximum"].asDouble();
			}

			if(json.isMember("exclusiveMinimum") && json["exclusiveMinimum"].isBool())
			{
				exclMin_ = json["exclusiveMinimum"].asBool();
			}

			if(json.isMember("exclusiveMaximum") && json["exclusiveMaximum"].isBool())
			{
				exclMax_ = json["exclusiveMaximum"].asBool();
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

			if(multSet_ && fmod(json.asDouble(), multipleOf_) != 0)
				PushError(path + ": Value must be a multiple of " + std::to_string(multipleOf_));

			if(minSet_)
			{
				if(exclMin_ && json.asDouble() <= min_)
					PushError(path + ": Value must be greater than " + std::to_string(min_));
				else if(json.asDouble() < min_)
					PushError(path + ": Value cannot be less than " + std::to_string(min_));
			}

			if(maxSet_)
			{
				if(exclMax_ && json.asDouble() >= max_)
					PushError(path + ": Value must be less than " + std::to_string(max_));
				else if(json.asDouble() > max_)
					PushError(path + ": Value cannot be greater than " + std::to_string(max_));
			}
		}
	};
}
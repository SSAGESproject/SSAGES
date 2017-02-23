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

#include <regex>
#include "Requirement.h"

namespace Json
{
	//! Requirements on strings.
	/*!
	 * \ingroup Json
	 */
	class StringRequirement : public Requirement 
	{
	private: 
		bool minSet_; //!< If \c True, minimum length requirement is active.
		bool maxSet_; //!< If \c True, maximum length requirement is active.
		bool rgxSet_; //!< If \c True, string has to match regular expression.
		size_t minLength_; //!< Minimum string length;
		size_t maxLength_; //!< Maximum string length;
		std::regex regex_; //!< Regular expression to match string to.
		std::string expr_; //!< Expression.
		std::string path_; //!< Path for JSON path specification.
		std::vector<std::string> enum_; //!< Enum values.

	public:
		//! Constructor.
		StringRequirement() : 
		minSet_(false), maxSet_(false), rgxSet_(false),
		minLength_(0), maxLength_(0), regex_(), expr_(), path_(), enum_(0)
		{}
		
		//! Reset Requirement.
		virtual void Reset() override
		{
			minSet_ = false;
			maxSet_ = false;
			rgxSet_ = false;
			minLength_ = 0;
			maxLength_ = 0;
			regex_ = "";
			expr_ = "";
			path_  = "";
			enum_.clear();
			ClearErrors();
			ClearNotices();
		}

		//! Parse JSON value to generate Requirement.
		/*!
		 * \param json JSON input value.
		 * \param path Path for JSON path specification.
		 */
		virtual void Parse(Value json, const std::string& path) override
		{
			Reset();
			
			path_ = path;
			if(json.isMember("minLength") && json["minLength"].isUInt())
			{
				minSet_ = true;
				minLength_ = json["minLength"].asUInt();
			}
			
			if(json.isMember("maxLength") && json["maxLength"].isUInt())
			{
				maxSet_ = true;
				maxLength_ = json["maxLength"].asUInt();
			
			}

			if(json.isMember("pattern") && json["pattern"].isString())
			{
				rgxSet_ = true;
				expr_ = json["pattern"].asString();
				regex_ = std::regex(expr_, std::regex::ECMAScript);
			}

			if(json.isMember("enum") && json["enum"].isArray())
			{
				for(const auto& val : json["enum"])
					enum_.push_back(val.asString());
			}
		}

		//! Validate string value.
		/*!
		 * \param json JSON value to be validated.
		 * \param path Path for JSON path specification.
		 *
		 * This function tests if the JSON value is of type string and if the
		 * string meets the requirements loaded via StringRequirement::Parse().
		 * If the validation fails, one or more errors are appended to the list
		 * of errors.
		 */
		virtual void Validate(const Value& json, const std::string& path) override
		{
			if(!json.isString())
			{
				PushError(path + ": Must be of type \"string\"");
				return;
			}
			
			if(minSet_ && json.asString().length() < minLength_)
				PushError(path + ": Length must be greater than " + std::to_string(minLength_));
			
			if(maxSet_ && json.asString().length() > maxLength_)
				PushError(path + ": Length must be less than " + std::to_string(minLength_));

			if(rgxSet_ && !std::regex_match(json.asString(), regex_))
				PushError(path + ": String must match regular expression \"" + expr_ + "\"");

			if(enum_.size() && std::find(enum_.begin(),enum_.end(), json.asString()) == enum_.end())
				PushError(path + ": String is not a valid entry");
		}
	};
}
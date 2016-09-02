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
		bool _minSet; //!< If \c True, minimum length requirement is active.
		bool _maxSet; //!< If \c True, maximum length requirement is active.
		bool _rgxSet; //!< If \c True, string has to match regular expression.
		size_t _minLength; //!< Minimum string length;
		size_t _maxLength; //!< Maximum string length;
		std::regex _regex; //!< Regular expression to match string to.
		std::string _expr; //!< Expression.
		std::string _path; //!< Path for JSON path specification.
		std::vector<std::string> _enum; //!< Enum values.

	public:
		//! Constructor.
		StringRequirement() : 
		_minSet(false), _maxSet(false), _rgxSet(false),
		_minLength(0), _maxLength(0), _regex(), _expr(), _path(), _enum(0)
		{}
		
		//! Reset Requirement.
		virtual void Reset() override
		{
			_minSet = false;
			_maxSet = false;
			_rgxSet = false;
			_minLength = 0;
			_maxLength = 0;
			_regex = "";
			_expr = "";
			_path  = "";
			_enum.clear();
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
			
			_path = path;
			if(json.isMember("minLength") && json["minLength"].isUInt())
			{
				_minSet = true;
				_minLength = json["minLength"].asUInt();
			}
			
			if(json.isMember("maxLength") && json["maxLength"].isUInt())
			{
				_maxSet = true;
				_maxLength = json["maxLength"].asUInt();
			
			}

			if(json.isMember("pattern") && json["pattern"].isString())
			{
				_rgxSet = true;
				_expr = json["pattern"].asString();
				_regex = std::regex(_expr, std::regex::ECMAScript);
			}

			if(json.isMember("enum") && json["enum"].isArray())
			{
				for(const auto& val : json["enum"])
					_enum.push_back(val.asString());
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
			
			if(_minSet && json.asString().length() < _minLength)
				PushError(path + ": Length must be greater than " + std::to_string(_minLength));
			
			if(_maxSet && json.asString().length() > _maxLength)
				PushError(path + ": Length must be less than " + std::to_string(_minLength));

			if(_rgxSet && !std::regex_match(json.asString(), _regex))
				PushError(path + ": String must match regular expression \"" + _expr + "\"");

			if(_enum.size() && std::find(_enum.begin(),_enum.end(), json.asString()) == _enum.end())
				PushError(path + ": String is not a valid entry");
		}
	};
}
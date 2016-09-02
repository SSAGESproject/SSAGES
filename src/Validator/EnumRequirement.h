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
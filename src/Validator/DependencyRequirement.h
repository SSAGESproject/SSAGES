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
	//! Requires dependencies to be met.
	/*!
	 * \ingroup Json
	 */
	class DependencyRequirement : public Requirement
	{
	private:
		//! Set of dependencies.
		std::map<std::string, std::vector<std::string>> _deps;

	public:
		//! Constructor
		DependencyRequirement() : _deps() {}

		//! Reset Requirement.
		/*!
		 * Clear set of dependencies.
		 */
		virtual void Reset() override
		{
			_deps.clear();
			ClearErrors();
			ClearNotices();
		}

		//! Parse JSON input value.
		/*!
		 * \param json JSON input value.
		 */
		virtual void Parse(Value json, const std::string&) override
		{
			Reset();

			auto names = json.getMemberNames();
			int i = 0;
			for(auto& value : json)
			{
				if(value.isArray())
				{
					_deps[names[i]] = {};
					for(auto& dep : value)
					{
						_deps[names[i]].push_back(dep.asString());
					}
				}
				++i;
			}
		}

		//! Validate that Requirement is met.
		/*!
		 * \param json JSON value to be validated.
		 * \param path Path for JSON path specification.
		 */
		virtual void Validate(const Value& json, const std::string& path) override
		{
			if(!json.isObject())
			{
				PushError(path + ": Dependency specified for non-object");
				return;
			}

			for(auto& dep : _deps)
			{
				auto& name = dep.first;
				if(json.isMember(name))
					for(const auto& val : dep.second)
						if(!json.isMember(val))
							PushError(path + ": \"" + name + "\" depends on \"" + val + "\"");
			}
		}
		
		//! Destructor.
		~DependencyRequirement() {}	
	};
}
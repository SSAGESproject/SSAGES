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
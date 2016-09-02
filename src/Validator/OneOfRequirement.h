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
#include "RequirementLoader.h"

namespace Json
{
	//! Requires exactly one of a list of Requirements to hold.
	/*!
	 * \ingroup Json
	 */
	class OneOfRequirement : public Requirement
	{
	private:
		RequireList _reqs; //!< List of Requirements.

	public:
		//! Constructor.
		OneOfRequirement() : _reqs(0) {}

		//! Clear list of error messages for all Requirements.
		virtual void ClearErrors() override
		{
			for(auto& r : _reqs)
				r->ClearErrors();

			Requirement::ClearErrors();
		}

		//! Clear list of notices for all Requirements.
		virtual void ClearNotices() override
		{
			for(auto& r : _reqs)
				r->ClearNotices();

			Requirement::ClearNotices();
		} 

		//! Reset requirement.
		virtual void Reset() override
		{
			ClearErrors();
			ClearNotices();
			_reqs.clear();
		}

		//! Parse JSON value to generate Requirement.
		/*!
		 * \param json JSON input value.
		 * \param path Path for JSON path specification.
		 */
		virtual void Parse(Value json, const std::string& path) override
		{
			Reset();
			RequirementLoader loader;

			auto& head = json.isMember("oneOf") ? json["oneOf"] : json;

			for(auto& val : head)
				if(auto req = loader.LoadRequirement(val))
				{
					_reqs.push_back(std::move(req));
					_reqs.back()->Parse(val, path);
				}

		}

		//! Validate that only one of the Requirements hold.
		/*!
		 * \param json JSON value to be validated.
		 * \param path Path for JSON path specification.
		 */
		virtual void Validate(const Value& json, const std::string& path) override
		{
			int validated = 0;
			for(auto& r : _reqs)
			{
				r->Validate(json, path);
				if(!r->HasErrors())
					++validated;
			}

			if(validated > 1)
				PushError(path + ": Input must validate against only one schema");
			else if(validated == 0)
				for(auto& r : _reqs)
				{
					if(r->HasErrors())
						for(const auto& error : r->GetErrors())
							PushError(error);
			
					if(r->HasNotices())
						for(const auto& notice : r->GetNotices())
							PushNotice(notice);
				}
		}
	};
}
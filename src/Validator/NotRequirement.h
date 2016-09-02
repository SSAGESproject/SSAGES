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
	//! Requires a given Requirement to fail.
	/*!
	 * \ingroup Json
	 */
	class NotRequirement : public Requirement
	{
	private:
		std::unique_ptr<Requirement> _req; //!< Requirement to negate.

	public:
		//! Constructor.
		NotRequirement() : _req(nullptr) {}

		//! Clear list of error messages.
		virtual void ClearErrors() override
		{
			if(_req)
				_req->ClearErrors();

			Requirement::ClearErrors();
		}

		//! Clear list of notices.
		virtual void ClearNotices() override
		{
			if(_req != nullptr)
				_req->ClearNotices();
		
			Requirement::ClearNotices();
		} 

		//! Reset Requirement.
		virtual void Reset() override
		{
			ClearErrors();
			ClearNotices();
			
			_req.reset();		
		}

		//! Parse JSON value and generate Requirement to be negated.
		/*!
		 * \param json JSON input value.
		 * \param path Path for JSON path specification.
		 */
		virtual void Parse(Value json, const std::string& path) override
		{
			Reset();
			RequirementLoader loader;

			auto& head = json.isMember("not") ? json["not"] : json;
			if((_req = loader.LoadRequirement(head)))
				_req->Parse(head, path);
		}

		//! Validate that JSON value fails the given Requirement.
		/*!
		 * \param json JSON value to validate.
		 * \param path Path for JSON path specification.
		 *
		 * Calls validate on the Requirement to be negated and checks if it has
		 * an error. If it has not, it passed the validation and an error is
		 * added to this Requirement.
		 */
		virtual void Validate(const Value& json, const std::string& path) override
		{
			_req->Validate(json, path);
			if(!_req->HasErrors())
				PushError(path + ": Value must not validate against requirement.");
		}
	};
}
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

#include <vector>
#include <map>
#include <memory>
#include "json/json.h"

namespace Json
{
	//! Requirements on input files
	/*!
	 * \ingroup Json
	 */
	class Requirement
	{
	private:
		std::vector<std::string> _errors; //!< List of error messages.
		std::vector<std::string> _notices; //!< List of messages.

	protected:
		//! Add error to list of error messages.
		/*!
		 * \param error Error message.
		 *
		 * This function adds an error message to the list of error messages.
		 * The list of error messages can be retrieved by Requirement::GetErrors()
		 */
		void PushError(const std::string& error) { _errors.push_back(error); }

		//! Add message to list of notices.
		/*!
		 * \param notice Message string.
		 *
		 * This function adds a new message to the list of messages. The list
		 * can be retrieved using Requirement::GetNotices().
		 */
		void PushNotice(const std::string& notice) { _notices.push_back(notice); }
	
	public:
		//! Parse JSON value.
		/*!
		 * \param json JSON value with input information.
		 * \param path Path for JSON path specification.
		 */
		virtual void Parse(Value json, const std::string& path) = 0;

		//! Validate that JSON value meets requirements.
		/*!
		 * \param json JSON value to validate.
		 * \param path Path for JSON path specification.
		 */
		virtual void Validate(const Value& json, const std::string& path) = 0;

		//! Reset validator
		virtual void Reset() = 0;

		//! Check if errors have occured.
		/*!
		 * \returns \c True if the list of errors is not empty.
		 */
		bool HasErrors() { return _errors.size() != 0; };

		//! Get list of error messages.
		/*!
		 * \return List of error messages.
		 */
		std::vector<std::string> GetErrors() { return _errors; };

		//! Clear list of error messages.
		virtual void ClearErrors() { _errors.clear(); }

		//! Check if notices have been queued.
		/*!
		 * \return \c True if list of notices is not empty.
		 */
		virtual bool HasNotices() {return _notices.size() != 0; };

		//! Get list of notices.
		/*!
		 * \return List of notice messages.
		 */
		std::vector<std::string> GetNotices() {return _notices; };

		//! Clear list of notice messages.
		virtual void ClearNotices() { _notices.clear(); }

		//! Destructor
		virtual ~Requirement() {}
	};

	//! List of Requirements
	using RequireList = std::vector<std::unique_ptr<Requirement>>;
}
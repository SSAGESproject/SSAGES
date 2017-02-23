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
		std::vector<std::string> errors_; //!< List of error messages.
		std::vector<std::string> notices_; //!< List of messages.

	protected:
		//! Add error to list of error messages.
		/*!
		 * \param error Error message.
		 *
		 * This function adds an error message to the list of error messages.
		 * The list of error messages can be retrieved by Requirement::GetErrors()
		 */
		void PushError(const std::string& error) { errors_.push_back(error); }

		//! Add message to list of notices.
		/*!
		 * \param notice Message string.
		 *
		 * This function adds a new message to the list of messages. The list
		 * can be retrieved using Requirement::GetNotices().
		 */
		void PushNotice(const std::string& notice) { notices_.push_back(notice); }
	
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
		bool HasErrors() { return errors_.size() != 0; };

		//! Get list of error messages.
		/*!
		 * \return List of error messages.
		 */
		std::vector<std::string> GetErrors() { return errors_; };

		//! Clear list of error messages.
		virtual void ClearErrors() { errors_.clear(); }

		//! Check if notices have been queued.
		/*!
		 * \return \c True if list of notices is not empty.
		 */
		virtual bool HasNotices() {return notices_.size() != 0; };

		//! Get list of notices.
		/*!
		 * \return List of notice messages.
		 */
		std::vector<std::string> GetNotices() {return notices_; };

		//! Clear list of notice messages.
		virtual void ClearNotices() { notices_.clear(); }

		//! Destructor
		virtual ~Requirement() {}
	};

	//! List of Requirements
	using RequireList = std::vector<std::unique_ptr<Requirement>>;
}
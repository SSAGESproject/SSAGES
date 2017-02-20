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
		std::unique_ptr<Requirement> req_; //!< Requirement to negate.

	public:
		//! Constructor.
		NotRequirement() : req_(nullptr) {}

		//! Clear list of error messages.
		virtual void ClearErrors() override
		{
			if(req_)
				req_->ClearErrors();

			Requirement::ClearErrors();
		}

		//! Clear list of notices.
		virtual void ClearNotices() override
		{
			if(req_ != nullptr)
				req_->ClearNotices();
		
			Requirement::ClearNotices();
		} 

		//! Reset Requirement.
		virtual void Reset() override
		{
			ClearErrors();
			ClearNotices();
			
			req_.reset();		
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
			if((req_ = loader.LoadRequirement(head)))
				req_->Parse(head, path);
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
			req_->Validate(json, path);
			if(!req_->HasErrors())
				PushError(path + ": Value must not validate against requirement.");
		}
	};
}
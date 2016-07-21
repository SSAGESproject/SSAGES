#pragma once 

#include "Requirement.h"
#include "RequirementLoader.h"

namespace Json
{
	//! Requires that at least one of a list of Requirements hold.
	/*!
	 * \ingroup Json
	 */
	class AnyOfRequirement : public Requirement
	{
	private:
		RequireList _reqs; //!< List of Requirements.

	public:
		//! Constructor.
		AnyOfRequirement() : _reqs(0) {}

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

		//! Reset all Requirements.
		virtual void Reset() override
		{
			ClearErrors();
			ClearNotices();
			_reqs.clear();
		}

		//! Parse JSON value and create list of Requirements.
		/*!
		 * \param json JSON input value.
		 * \param path Path for JSON path specification.
		 *
		 * This function parses the given JSON value and creates a new
		 * Requirement for each value in the "anyOf" branch. The new
		 * Requirements are appended to the list of Requirements.
		 */
		virtual void Parse(Value json, const std::string& path) override
		{
			Reset();
			RequirementLoader loader;

			auto& head = json.isMember("anyOf") ? json["anyOf"] : json;

			for(auto& val : head)
				if(auto req = loader.LoadRequirement(val))
				{
					_reqs.push_back(std::move(req));
					_reqs.back()->Parse(val, path);
				}

		}

		//! Validate that at least one Requirement holds.
		/*!
		 * \param json JSON value to validate.
		 * \param path Path for JSON path specification.
		 */
		virtual void Validate(const Value& json, const std::string& path) override
		{
			for(auto& r : _reqs)
			{
				r->Validate(json, path);
				if(!r->HasErrors())
					return;
			}

			// Collect errors.
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
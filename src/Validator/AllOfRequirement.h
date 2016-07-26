#pragma once 

#include "Requirement.h"
#include "RequirementLoader.h"
#include <memory>

namespace Json
{
	//! Requires that all of a list of Requirements hold.
	/*!
	 * \ingroup Json
	 */
	class AllOfRequirement : public Requirement
	{
	private:
		RequireList _reqs; //!< List of Requirements.

	public:
		//! Constructor
		AllOfRequirement() : _reqs(0) {}

		//! Clear lists of errors for all Requirements.
		virtual void ClearErrors() override
		{
			for(auto& r : _reqs)
				r->ClearErrors();

			Requirement::ClearErrors();
		}

		//! Clear lists of notices for all Requirements.
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

		//! Parse JSON value and add all loaded Requirements.
		/*!
		 * \param json JSON input value.
		 * \param path Path for JSON path specification.
		 *
		 * This function parses a JSON value and adds all Requirements that are
		 * in the "allOf" branch to the list of Requirements.
		 */
		virtual void Parse(Value json, const std::string& path) override
		{
			Reset();
			RequirementLoader loader;

			auto& head = json.isMember("allOf") ? json["allOf"] : json;

			for(auto& val : head)
				if(auto req = loader.LoadRequirement(val))
				{
					_reqs.push_back(std::move(req));
					_reqs.back()->Parse(val, path);
				}

		}

		//! Validate that all Requirements are met.
		/*!
		 * \param json JSON value to validate.
		 * \param path Path for JSON path specification.
		 */
		virtual void Validate(const Value& json, const std::string& path) override
		{
			for(auto& r : _reqs)
			{
				r->Validate(json, path);
		
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
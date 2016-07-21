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
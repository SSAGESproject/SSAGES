#pragma once 

#include "Requirement.h"

namespace Json
{
	//! Requires value to be of type Null.
	/*!
	 * \ingroup Json
	 */
	class NullRequirement : public Requirement
	{
	public:
		//! Reset this Requirement.
		virtual void Reset() {}

		//! Parse JSON value to set up this Requirement.
		virtual void Parse(Value, const std::string&) {}

		//! Validate that JSON value is null.
		/*!
		 * \param json JSON value to validate.
		 * \param path Path for JSON path specification.
		 */
		virtual void Validate(const Value& json, const std::string& path)
		{
			if(!json.isNull())
				PushError(path + ": Must be a null value");
		}
	};
}
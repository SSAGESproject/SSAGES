#pragma once 

#include "Requirement.h"

namespace Json
{
	//! Requires json value to be of type Bool.
	/*!
	 * \ingroup Json
	 */
	class BooleanRequirement : public Requirement
	{
	public:
		//! Reset Requirement.
		virtual void Reset() {}

		//! Parse JSON string
		virtual void Parse(Value, const std::string&) {}

		//! Validate that JSON value is Bool.
		/*!
		 * \param json JSON value to validate.
		 * \param path Path for JSON path specification.
		 */
		virtual void Validate(const Value& json, const std::string& path)
		{
			if(!json.isBool())
				PushError(path + ": Must be a boolean");
		}
	};
}
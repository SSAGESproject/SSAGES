#pragma once 

#include "Requirement.h"

namespace Json
{
	class BooleanRequirement : public Requirement
	{
	public:
		virtual void Reset() {}

		virtual void Parse(Value, const std::string&) {}

		virtual void Validate(const Value& json, const std::string& path)
		{
			if(!json.isBool())
				PushError(path + ": Must be a boolean");
		}
	};
}
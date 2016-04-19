#pragma once 

#include "json/json.h"
#include "Requirement.h"
#include <memory>

namespace Json
{
	class RequirementLoader
	{
	public:
		std::unique_ptr<Requirement> LoadRequirement(const Value& json);

		std::unique_ptr<Requirement> LoadExtended(const Value& json);
		
		std::unique_ptr<Requirement> LoadRequirement(const ValueType& type);
	};
}
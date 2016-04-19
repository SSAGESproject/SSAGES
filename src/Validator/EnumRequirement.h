#pragma once 

#include "Requirement.h"

namespace Json
{
	class EnumRequirement : public Requirement
	{
	private:
		std::vector<Value> _enum;
	
	public:
		EnumRequirement() : _enum(0) {}

		virtual void Reset() 
		{
			_enum.clear();
		}

		virtual void Parse(Value json, const std::string&) 
		{
			if(json.isArray())
				for(auto& val : json)
					_enum.push_back(val);
		}

		virtual void Validate(const Value& json, const std::string& path)
		{
			bool found = false;

			for(auto& val : _enum)
				if(json == val)
					found = true;

			if(!found)
				PushError(path  + ": Value is not a valid entry.");
		}
	};
}
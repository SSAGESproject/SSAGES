#pragma once 

#include "Requirement.h"

namespace Json
{
	class DependencyRequirement : public Requirement
	{
	private:
		std::map<std::string, std::vector<std::string>> _deps;

	public:
		DependencyRequirement() : _deps() {}

		virtual void Reset() override
		{
			_deps.clear();
			ClearErrors();
			ClearNotices();
		}

		virtual void Parse(Value json, const std::string&) override
		{
			Reset();

			auto names = json.getMemberNames();
			int i = 0;
			for(auto& value : json)
			{
				if(value.isArray())
				{
					_deps[names[i]] = {};
					for(auto& dep : value)
					{
						_deps[names[i]].push_back(dep.asString());
					}
				}
				++i;
			}
		}

		virtual void Validate(const Value& json, const std::string& path) override
		{
			if(!json.isObject())
			{
				PushError(path + ": Dependency specified for non-object");
				return;
			}

			for(auto& dep : _deps)
			{
				auto& name = dep.first;
				if(json.isMember(name))
					for(const auto& val : dep.second)
						if(!json.isMember(val))
							PushError(path + ": \"" + name + "\" depends on \"" + val + "\"");
			}
		}
		
		~DependencyRequirement() {}	
	};
}
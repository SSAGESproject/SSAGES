#pragma once 

#include <iostream>
#include <algorithm>
#include <list>
#include <map>
#include "Requirement.h"
#include "StringRequirement.h"
#include "IntegerRequirement.h"
#include "NumberRequirement.h"
#include "DependencyRequirement.h"
#include "RequirementLoader.h"

namespace Json
{
	class ObjectRequirement : public Requirement
	{
	private:
		std::map<std::string, std::unique_ptr<Requirement>> _properties;
		std::map<std::string, std::unique_ptr<Requirement>> _patternProps;
		RequireList _extended;
		std::unique_ptr<DependencyRequirement> _dependency;

		std::vector<std::string> _required;
		bool _moreProps, _setMin, _setMax;
		unsigned int _min, _max;

	public:
		ObjectRequirement() : 
		_properties(), _patternProps(), _extended(0), _dependency(nullptr), _required(),
		_moreProps(true), _setMin(false), _setMax(false), _min(0), _max(0)
		{}

		~ObjectRequirement()
		{
			_properties.clear();

			_patternProps.clear();

			_extended.clear();
		}

		virtual void ClearErrors() override
		{
			for(auto& c : _properties)
				c.second->ClearErrors();

			for(auto& c : _patternProps)
				c.second->ClearErrors();

			for(auto& c : _extended)
				c->ClearErrors();

			if(_dependency != nullptr)
				_dependency->ClearErrors();

			Requirement::ClearErrors();
		}

		virtual void ClearNotices() override
		{
			for(auto& c : _properties)
				c.second->ClearNotices();

			for(auto& c : _patternProps)
				c.second->ClearNotices();

			for(auto& c : _extended)
				c->ClearNotices();

			if(_dependency != nullptr)
				_dependency->ClearNotices();

			Requirement::ClearNotices();
		} 

		virtual void Reset() override
		{
			ClearErrors();
			ClearNotices();

			_properties.clear();
			_patternProps.clear();

			_moreProps = true;
			_setMin = _setMax = false;
			_min = _max = 0;
			_required.clear();
			_dependency.reset();
		}

		virtual void Parse(Value json, const std::string& path) override
		{
			Reset();
			RequirementLoader loader;

			// Additional properties.
			if(json.isMember("additionalProperties"))
			{
				if(json["additionalProperties"].isBool())
					_moreProps = json["additionalProperties"].asBool();
				else if(json["additionalProperties"].isObject())
					json["properties"]["additionalProperties"] = json["additionalProperties"];
			}

			// properties.
			if(json.isMember("properties") && json["properties"].isObject())
			{
				auto& props = json["properties"];
				auto names = props.getMemberNames();
				int i = 0;
				for(auto& prop : props)
				{
					if(auto property = loader.LoadRequirement(prop))
					{
						_properties[names[i]] = std::move(property);
						_properties[names[i]]->Parse(prop, path + "/" + names[i]);
					}

					++i;
				}
			}

			// Pattern properties. TODO: eliminate redundant code!!
			if(json.isMember("patternProperties") && json["patternProperties"].isObject())
			{
				auto& props = json["patternProperties"];
				auto names = props.getMemberNames();
				int i = 0;
				for(auto& prop : props)
				{
					if(prop.isObject())
					{
						if(auto property = loader.LoadRequirement(prop))
						{
							_patternProps[names[i]] = std::move(property);
							_patternProps[names[i]]->Parse(prop, path + "/" + names[i]);
						}
					}

					++i;
				}
			}

			// Required properties.
			if(json.isMember("required") && json["required"].isArray())
			{
				for(auto& requirement : json["required"])
					_required.push_back(requirement.asString());
			}

			// Min property count.
			if(json.isMember("minProperties") && json["minProperties"].isUInt())
			{
				_setMin = true;
				_min = json["minProperties"].asInt();
			}

			// Max property count.
			if(json.isMember("maxProperties") && json["maxProperties"].isUInt())
			{
				_setMax = true;
				_max = json["maxProperties"].asInt();
			}

			// Dependencies
			if(json.isMember("dependencies") && json["dependencies"].isObject())
			{
				_dependency = std::move(std::unique_ptr<DependencyRequirement>(new DependencyRequirement()));
				_dependency->Parse(json["dependencies"], path);
			}

			// Extended properties. 
			for(auto& prop : json)
			{
				if(auto req  = loader.LoadExtended(prop))
				{
					_extended.push_back(std::move(req));
					_extended.back()->Parse(prop, path);
				}
			}
		}

		virtual void Validate(const Value& json, const std::string& path) override
		{
			if(!json.isObject())
			{
				PushError(path + ": Must be of type \"object\"");
				return;
			}

			if(_setMin && json.size() < _min)
				PushError(path + ": Object must contain at least " + std::to_string(_min) + " properties");

			if(_setMax && json.size() > _max)
				PushError(path + ": Object must contain at most " + std::to_string(_max) + " properties");


			// Check dependencies. 
			if(_dependency != nullptr)
			{
				_dependency->Validate(json, path);
				if(_dependency->HasErrors())
					for(const auto& error : _dependency->GetErrors())
						PushError(error);
				if(_dependency->HasNotices())
					for(const auto& notice : _dependency->GetNotices())
						PushNotice(notice);
			}

			// Copy so we can pop items off the list.
			auto rprops = _required;

			auto names = json.getMemberNames();
			int i = 0;
			for(auto& prop : json)
			{
				Requirement* requirement = nullptr; 
				auto it = _properties.find(names[i]);
				if(it != _properties.end())
					requirement = it->second.get();					
				else if(_patternProps.size() != 0)
				{
					for(auto& pattern : _patternProps)
					{
						auto regex = std::regex(pattern.first, std::regex::ECMAScript);
						if(std::regex_search(names[i], regex))
							requirement = pattern.second.get();
					}
				}
				
				if(!requirement && _properties.find("additionalProperties") != _properties.end())
					requirement = _properties["additionalProperties"].get();
				
				if(requirement)
				{
					requirement->Validate(prop, path + "/" + names[i]);
					if(requirement->HasErrors())
						for(const auto& error : requirement->GetErrors())
							PushError(error);
					if(requirement->HasNotices())
						for(const auto& notice : requirement->GetNotices())
							PushNotice(notice);
				}
				else if(!_moreProps)
					PushError(path + ": Invalid property \"" + names[i] + "\" specified");

				rprops.erase(std::remove(rprops.begin(), rprops.end(),names[i]),rprops.end());
				++i;
			}

			if(_required.size() && rprops.size() != 0)
			{
				std::string msg = std::accumulate(rprops.begin(), rprops.end(), std::string(), 
    				[](const std::string& a, const std::string& b) -> std::string { 
        				return a + (a.length() > 0 ? ", " : "") + b; 
    			});
				PushError(path + ": Missing properties: " + msg);
			}

			// Validate extended.
			for(auto& requirement : _extended)
			{
				requirement->Validate(json, path);
				if(requirement->HasErrors())
						for(const auto& error : requirement->GetErrors())
							PushError(error);
				if(requirement->HasNotices())
					for(const auto& notice : requirement->GetNotices())
						PushNotice(notice);
			}
		}
	};
}
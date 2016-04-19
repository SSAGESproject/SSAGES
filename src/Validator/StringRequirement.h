#pragma once 

#include <regex>
#include "Requirement.h"

namespace Json
{
	class StringRequirement : public Requirement 
	{
	private: 
		bool _minSet, _maxSet, _rgxSet;
		size_t _minLength, _maxLength;
		std::regex _regex;
		std::string _expr;
		std::string _path;
		std::vector<std::string> _enum;

	public:
		StringRequirement() : 
		_minSet(false), _maxSet(false), _rgxSet(false),
		_minLength(0), _maxLength(0), _regex(), _expr(), _path(), _enum(0)
		{}
		
		virtual void Reset() override
		{
			_minSet = false;
			_maxSet = false;
			_rgxSet = false;
			_minLength = 0;
			_maxLength = 0;
			_regex = "";
			_expr = "";
			_path  = "";
			_enum.clear();
			ClearErrors();
			ClearNotices();
		}

		virtual void Parse(Value json, const std::string& path) override
		{
			Reset();
			
			_path = path;
			if(json.isMember("minLength") && json["minLength"].isUInt())
			{
				_minSet = true;
				_minLength = json["minLength"].asUInt();
			}
			
			if(json.isMember("maxLength") && json["maxLength"].isUInt())
			{
				_maxSet = true;
				_maxLength = json["maxLength"].asUInt();
			
			}

			if(json.isMember("pattern") && json["pattern"].isString())
			{
				_rgxSet = true;
				_expr = json["pattern"].asString();
				_regex = std::regex(_expr, std::regex::ECMAScript);
			}

			if(json.isMember("enum") && json["enum"].isArray())
			{
				for(const auto& val : json["enum"])
					_enum.push_back(val.asString());
			}
		}

		virtual void Validate(const Value& json, const std::string& path) override
		{
			if(!json.isString())
			{
				PushError(path + ": Must be of type \"string\"");
				return;
			}
			
			if(_minSet && json.asString().length() < _minLength)
				PushError(path + ": Length must be greater than " + std::to_string(_minLength));
			
			if(_maxSet && json.asString().length() > _maxLength)
				PushError(path + ": Length must be less than " + std::to_string(_minLength));

			if(_rgxSet && !std::regex_match(json.asString(), _regex))
				PushError(path + ": String must match regular expression \"" + _expr + "\"");

			if(_enum.size() && std::find(_enum.begin(),_enum.end(), json.asString()) == _enum.end())
				PushError(path + ": String is not a valid entry");
		}
	};
}
#pragma once 

#include <cmath>
#include "Requirement.h"

namespace Json
{
	class NumberRequirement : public Requirement
	{
	private:
		std::string _path;
		double _multipleOf, _min, _max;
		bool _multSet, _minSet, _maxSet, _exclMin, _exclMax;


	public:
		NumberRequirement() : 
		_path(), _multipleOf(0), _min(0), _max(0), _multSet(false), 
		_minSet(false), _maxSet(false), _exclMin(false), _exclMax(false)
		{}

		virtual void Reset() override
		{
			_multipleOf = 0;
			_minSet = _maxSet = false;
			_exclMin = _exclMax = false; 
			_min = _max = 0;
			_multSet = false;
			ClearErrors();
			ClearNotices();
		}

		virtual void Parse(Value json, const std::string& path) override
		{
			Reset();
			
			_path = path;
			if(json.isMember("multipleOf") && json["multipleOf"].isNumeric())
			{
				_multSet = true;
				_multipleOf = json["multipleOf"].asDouble();
			}

			if(json.isMember("minimum") && json["minimum"].isNumeric())
			{
				_minSet = true;
				_min = json["minimum"].asDouble();
			}

			if(json.isMember("maximum") && json["maximum"].isNumeric())
			{
				_maxSet = true;
				_max = json["maximum"].asDouble();
			}

			if(json.isMember("exclusiveMinimum") && json["exclusiveMinimum"].isBool())
			{
				_exclMin = json["exclusiveMinimum"].asBool();
			}

			if(json.isMember("exclusiveMaximum") && json["exclusiveMaximum"].isBool())
			{
				_exclMax = json["exclusiveMaximum"].asBool();
			}
		}

		virtual void Validate(const Value& json, const std::string& path) override
		{
			if(!json.isNumeric())
			{
				PushError(path + ": Must be of type \"number\".");
				return;
			}

			if(_multSet && fmod(json.asDouble(), _multipleOf) != 0)
				PushError(path + ": Value must be a multiple of " + std::to_string(_multipleOf));

			if(_minSet)
			{
				if(_exclMin && json.asDouble() <= _min)
					PushError(path + ": Value must be greater than " + std::to_string(_min));
				else if(json.asDouble() < _min)
					PushError(path + ": Value cannot be less than " + std::to_string(_min));
			}

			if(_maxSet)
			{
				if(_exclMax && json.asDouble() >= _max)
					PushError(path + ": Value must be less than " + std::to_string(_max));
				else if(json.asDouble() > _max)
					PushError(path + ": Value cannot be greater than " + std::to_string(_max));
			}
		}
	};
}
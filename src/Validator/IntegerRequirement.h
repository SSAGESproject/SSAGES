#pragma once 

#include "Requirement.h"

namespace Json
{
	class IntegerRequirement : public Requirement
	{
	private:
		std::string _path;
		int _multipleOf, _min, _max;
		bool _multSet, _minSet, _maxSet, _exclMin, _exclMax;

	public:
		IntegerRequirement() : 
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
			if(json.isMember("multipleOf") && json["multipleOf"].isInt())
			{
				_multSet = true;
				_multipleOf = json["multipleOf"].asInt();
			}

			if(json.isMember("minimum") && json["minimum"].isInt())
			{
				_minSet = true;
				_min = json["minimum"].asInt();
			}

			if(json.isMember("maximum") && json["maximum"].isInt())
			{
				_maxSet = true;
				_max = json["maximum"].asInt();
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
			if(!json.isInt())
			{
				PushError(path + ": Must be of type \"integer\".");
				return;
			}

			if(_multSet && (json.asInt() % _multipleOf != 0))
				PushError(path + ": Value must be a multiple of " + std::to_string(_multipleOf));

			if(_minSet)
			{
				if(_exclMin && json.asInt() <= _min)
					PushError(path + ": Value must be greater than " + std::to_string(_min));
				else if(json.asInt() < _min)
					PushError(path + ": Value cannot be less than " + std::to_string(_min));
			}

			if(_maxSet)
			{
				if(_exclMax && json.asInt() >= _max)
					PushError(path + ": Value must be less than " + std::to_string(_max));
				else if(json.asInt() > _max)
					PushError(path + ": Value cannot be greater than " + std::to_string(_max));
			}
		}
	};
}
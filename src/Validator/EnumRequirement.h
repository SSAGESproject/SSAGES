#pragma once 

#include "Requirement.h"

namespace Json
{
	//! Requires entry to be member of an enum.
	/*!
	 * \ingroup Json
	 */
	class EnumRequirement : public Requirement
	{
	private:
		std::vector<Value> _enum; //!< Enum value.
	
	public:
		//! Constructor
		EnumRequirement() : _enum(0) {}

		//! Clear enum value.
		virtual void Reset() 
		{
			_enum.clear();
		}

		//! Parse JSON input value to generate enum.
		/*!
		 * \param json JSON input value.
		 */
		virtual void Parse(Value json, const std::string&) 
		{
			if(json.isArray())
				for(auto& val : json)
					_enum.push_back(val);
		}

		//! Validate that JSON value is member of the parsed enum.
		/*!
		 * \param json JSON value to be validated.
		 * \param path Path for JSON path specification.
		 */
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
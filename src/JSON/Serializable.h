#pragma once 

#include "json/json.h"

namespace SSAGES
{
	class Serializable
	{

	public:
		// Serialize the class.
		virtual void Serialize(Json::Value& json) const = 0;
	};
}
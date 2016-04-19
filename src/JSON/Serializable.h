#pragma once 

#include "json/json.h"

namespace SAPHRON
{
	class Serializable
	{

	public:
		// Serialize the class.
		virtual void Serialize(Json::Value& json) const = 0;
	};
}
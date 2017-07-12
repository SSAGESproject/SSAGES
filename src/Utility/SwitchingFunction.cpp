/**
 * This file is part of
 * SSAGES - Suite for Advanced Generalized Ensemble Simulations
 *
 * Copyright 2017 Hythem Sidky <hsidky@nd.edu>
 *
 * SSAGES is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * SSAGES is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with SSAGES.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "SwitchingFunction.h"
#include <stdexcept>

namespace SSAGES
{
	SwitchingFunction* SwitchingFunction::Build(const Json::Value& json)
	{
		auto type = json.get("type", "rational").asString(); 
		if(type == "rational")
			return RationalSF::Build(json);
		else
			throw std::invalid_argument("Invalid switching function type \"" + type + "\".");
	}
}
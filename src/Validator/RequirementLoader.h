/**
 * This file is part of
 * SSAGES - Suite for Advanced Generalized Ensemble Simulations
 *
 * Copyright 2016 Ben Sikora <bsikora906@gmail.com>
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
#pragma once 

#include "json/json.h"
#include "Requirement.h"
#include <memory>

namespace Json
{
    //! Helper class to load Requirement
    /*!
     * \ingroup Json
     */
	class RequirementLoader
	{
	public:
        //! Load specific requirement.
        /*!
         * \param json JSON value specifying the type of Requirement.
         * \return Unique pointer to the Requirement.
         */
		std::unique_ptr<Requirement> LoadRequirement(const Value& json);

        //! Extended Requirement loader.
        /*!
         * \param json JSON value specifying the type of Requirement.
         * \return Unique pointer to the Requirement.
         */
		std::unique_ptr<Requirement> LoadExtended(const Value& json);
		
        //! Load Requirement based on Value Type.
        /*!
         * \param type Value Type specifying the Requirement.
         * \return Unique pointer to the Requirement.
         *
         * Return a Requirement based on the given Value Type.
         */
		std::unique_ptr<Requirement> LoadRequirement(const ValueType& type);
	};
}
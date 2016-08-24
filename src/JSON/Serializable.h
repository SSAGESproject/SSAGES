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

namespace SSAGES
{
    //! Base class for all serializable objects.
    /*!
     * \ingroup Json
     */
    class Serializable
    {
        public:

        //! Serialize the class.
        /*!
         * \param json JSON value to which the information is stored
         *
         * Store all information pertaining to this class in the given JSON
         * value. The information can then be used to restart or analyze the
         * simulation
         */
        virtual void Serialize(Json::Value& json) const = 0;

        //! Destructor
        virtual ~Serializable() { }
    };
}
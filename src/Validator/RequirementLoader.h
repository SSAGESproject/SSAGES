/**
 * This file has been obtained from
 * SAPHRON - Statistical Applied PHysics through Random On-the-fly Numerics
 * https://github.com/hsidky/SAPHRON
 *
 * Copyright 2016 Hythem Sidky
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE. 
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
/**
 * This file has been adapted from
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
#include <mpi.h>

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

    //! Base class for JSON buildable objects. 
    /*!
     * \ingroup JSON 
     */
    template<class T> 
    class Buildable
    {
        //! Build an object from JSON node. 
		/*! 
		 * \param JSON Value containing all input information.
		 * \param world MPI global communicator. 
         * \param comm MPI local communicator. 
		 * \param path Path for JSON specification. 
		 * \return Pointer to built object.
		 *
		 * \note This function builds an object from a JSON node. It will generally 
		 *       throw an exception if an error occurs, but may also return nullptr 
		 *       if an unknown error occurs. 
		 *       Object lifetime is the caller's responsibility! 
		 */
        static T* Build(const Json::Value& json, 
                        const MPI_Comm& world, 
                        const MPI_Comm& comm, 
                        const std::string& path)
        {
            return T::Construct(json, world, comm, path);
        }
    };

}
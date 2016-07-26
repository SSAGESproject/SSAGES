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
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
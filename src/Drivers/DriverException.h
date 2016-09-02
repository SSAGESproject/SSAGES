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

#include <vector>
#include <iostream>
#include <exception>
#include <stdexcept>
#include <sstream>

namespace SSAGES
{
	//! Exception to be thrown when building the Driver fails.
	class BuildException : public std::runtime_error
	{
	private:
		std::vector<std::string> _errors; //!< Error message.
	public:
		//! Constructor
		/*!
		 * \param errors Error message
		 */
		BuildException(std::vector<std::string> errors) : 
			std::runtime_error("Object build error"), _errors(errors)
		{
		}

		//! Get generic error message
		/*!
		 * \return Generic error message.
		 */
		virtual const char* what() const throw()
		{
			std::ostringstream msg(""); 

			msg << runtime_error::what() << ": "
				<< "See errors for details.";

			return msg.str().c_str();
		}

		//! Get specific error message
		/*!
		 * \return Specific error message.
		 */
		std::vector<std::string> GetErrors()
		{
			return _errors;
		}
	};
}
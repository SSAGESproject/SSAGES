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
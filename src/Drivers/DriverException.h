#pragma once 

#include <vector>
#include <iostream>
#include <exception>
#include <stdexcept>
#include <sstream>

namespace SSAGES
{
	class BuildException : public std::runtime_error
	{
	private:
		std::vector<std::string> _errors;
	public:
		BuildException(std::vector<std::string> errors) : 
			std::runtime_error("Object build error"), _errors(errors)
		{
		}

		virtual const char* what() const throw()
		{
			std::ostringstream msg(""); 

			msg << runtime_error::what() << ": "
				<< "See errors for details.";

			return msg.str().c_str();
		}

		std::vector<std::string> GetErrors()
		{
			return _errors;
		}
	};
}
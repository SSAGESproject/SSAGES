/**
 * This file has been obtained from
 * SAPHRON - Statistical Applied PHysics through Random On-the-fly Numerics
 * https://github.com/hsidky/SAPHRON
 *
 * Copyright 2017 Hythem Sidky
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

#include <vector>
#include <iostream>
#include <exception>
#include <stdexcept>
#include <sstream>
#include <iomanip>
#include <unistd.h>

namespace SSAGES
{
	//! Exception to be thrown when building the Driver fails.
	class BuildException : public std::runtime_error
	{
	private:
		std::vector<std::string> errors_; //!< Error message.
	public:
		//! Constructor
		/*!
		 * \param errors Error message
		 */
		BuildException(std::vector<std::string> errors) : 
			std::runtime_error("Object build error"), errors_(errors)
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

			return strdup(msg.str().c_str());
		}

		//! Get specific error message
		/*!
		 * \return Specific error message.
		 */
		std::vector<std::string> GetErrors()
		{
			return errors_;
		}
	};

	//! Print out Notice in bold text
	/*!
	 * \param notice Text to be written to console.
	 * \param msgw Width of the message.
	 */
	inline void PrintBoldNotice(const std::string& notice, int msgw)
	{
		if(isatty(fileno(stdout)))
		{
			std::cout << std::setw(msgw + 8) << std::left << "\033[1m" + notice + "\033[0m";
		}
		else
		{
			std::cout << std::setw(msgw + 8) << std::left << notice;
		}
	}

	//! Print a list of errors
	/*!
	 * \param msgs List of error messages.
	 * \param notw Width of the messages.
	 * \return -1 always.
	 *
	 * Write a list of error messages to the console.
	 */
	inline int DumpErrorsToConsole(const std::vector<std::string>& msgs, int notw)
	{
		if(isatty(fileno(stdout)))
		{
			std::cout << std::setw(notw) << std::right << "\033[1;31mError(s)! See below.\033[0m\n";
		}
		else
		{
			std::cout << std::setw(notw) << std::right << "Error(s)! See below.\n";
		}
		for(auto& msg : msgs)
				std::cout << " * " << msg << "\n";
		return -1;
	}

	//! Print a list of notices
	/*!
	 * \param msgs List of messages.
	 * \param prefix Prefix to prepend to each message.
	 * \param notw Width of the messages.
	 */
	inline void DumpNoticesToConsole(const std::vector<std::string>& msgs, std::string prefix, int notw)
	{
		if(isatty(fileno(stdout)))
		{
			std::cout << std::setw(notw) << std::right << "\033[32mOK!\033[0m\n";
		}
		else
		{
			std::cout << std::setw(notw) << std::right << "OK!\n";
		}
		if(msgs.size() == 0)
			return;
		
		for(auto& msg : msgs)
			std::cout << prefix << " * " << msg << "\n";
	}
}
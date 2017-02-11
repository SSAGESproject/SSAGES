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
#include <iomanip>
#include <boost/mpi.hpp>

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

			return msg.str().c_str();
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
	 * \param world MPI global communicator.
	 */
	inline void PrintBoldNotice(const std::string& notice, int msgw, const boost::mpi::communicator& world)
	{

		if(world.rank() == 0)
			std::cout << std::setw(msgw + 8) << std::left << "\033[1m" + notice + "\033[0m";
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
		std::cout << std::setw(notw) << std::right << "\033[1;31mError(s)! See below.\033[0m\n";
		for(auto& msg : msgs)
				std::cout << " * " << msg << "\n";
		return -1;
	}

	//! Print a list of notices
	/*!
	 * \param msgs List of messages.
	 * \param prefix Prefix to prepend to each message.
	 * \param notw Width of the messages.
	 * \param world MPI global communicator.
	 */
	inline void DumpNoticesToConsole(const std::vector<std::string>& msgs, std::string prefix, int notw, const boost::mpi::communicator& world)
	{
		if(world.rank() == 0)
		{
			std::cout << std::setw(notw) << std::right << "\033[32mOK!\033[0m\n";
			if(msgs.size() == 0)
				return;
			
			for(auto& msg : msgs)
				std::cout << prefix << " * " << msg << "\n";
		}
	}
}
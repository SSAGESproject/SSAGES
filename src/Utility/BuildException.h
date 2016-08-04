#pragma once

#include <vector>
#include <boost/mpi.hpp>
#include <iomanip>

namespace mpi = boost::mpi;
namespace SSAGES
{
	//! Exception to be thrown when setting up the Simulation, Driver, CV, etc.
	class BuildException : public std::runtime_error
	{
	private:
		std::vector<std::string> _errors; //!< Error string
	public:
		//! Construtor
		/*!
		 * \param errors Error message.
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

	//! Print out Notice in bold text
	/*!
	 * \param notice Text to be written to console.
	 * \param msgw Width of the message.
	 * \param world MPI global communicator.
	 */
	inline void PrintBoldNotice(const std::string& notice, int msgw, const mpi::communicator& world)
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
	inline void DumpNoticesToConsole(const std::vector<std::string>& msgs, std::string prefix, int notw, const mpi::communicator& world)
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

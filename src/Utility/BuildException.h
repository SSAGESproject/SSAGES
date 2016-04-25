#pragma once

#include <vector>
#include <boost/mpi.hpp>
#include <iomanip>

namespace mpi = boost::mpi;
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

	void PrintBoldNotice(const std::string& notice, int msgw, const mpi::communicator& world)
	{

		if(world.rank() == 0)
			std::cout << std::setw(msgw + 8) << std::left << "\033[1m" + notice + "\033[0m";
	}

	int DumpErrorsToConsole(const std::vector<std::string>& msgs, int notw)
	{
		std::cout << std::setw(notw) << std::right << "\033[1;31mError(s)! See below.\033[0m\n";
		for(auto& msg : msgs)
				std::cout << " * " << msg << "\n";
		return -1;
	}

	void DumpNoticesToConsole(const std::vector<std::string>& msgs, std::string prefix, int notw, const mpi::communicator& world)
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

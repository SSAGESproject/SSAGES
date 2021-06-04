/**
 * This file is part of
 * SSAGES - Software Suite for Advanced General Ensemble Simulations
 *
 * Copyright 2017 Hythem Sidky <hsidky@nd.edu>
 *                Ben Sikora <bsikora906@gmail.com>
 *                Julian Helfferich <julian.helfferich@gmail.com>
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
#include <getopt.h>
#include <string>
#include <sstream>

#include "config.h"
#include <mxx/comm.hpp>
#include <mxx/env.hpp>
#include "Driver.h"
#include "JSON/JSONLoader.h"
#include "Drivers/DriverException.h"

using namespace SSAGES;
using namespace Json;

int main(int argc, char* argv[])
{
	// Let MPI collect command line arguments pertaining to MPI.
	mxx::env env(argc, argv);
	mxx::env::set_exception_on_error();
	auto world = mxx::comm();

	std::stringstream helpStream;
	helpStream
		<< "Welcome to SSAGES - Software Suite for Advanced General Ensemble Simulations\n"
		<< "\n"
		<< argv[0] << " [-h | --help] [ --version ] InputFile.json\n"
		<< "\n"
		<< "Required parameters:\n"
		<< "InputFile.json - Input file specifying the general ensemble method and\n"
		<< "                 all required input parameters in the JSON format. For more\n"
		<< "                 information on the input file, confer the user manual.\n"
		<< "\n"
		<< "Optional parameters:\n"
		<< " -h | --help - Print this help message\n"
		<< " --version   - Print SSAGES and engine version and quit\n";

	std::stringstream versionStream;
	versionStream
		<< "SSAGES version " << SSAGES_VERSION << "\n"
		<< "Engine: " << MD_ENGINE_VERSION << "\n";

	// Define short options, e.g. -h
	const char *OPT_STRING = "h!";

	// Define long options, e.g. --help
	static struct option longOpts[] = {
		{ "help", no_argument, nullptr, 'h' },
		{ "version", no_argument, nullptr, '!' },
	};

	// Parse command line arguments
	int optionindex = 1;
	int opt = getopt_long(argc, argv, OPT_STRING, longOpts, &optionindex);

	while (opt != -1) {
		switch (opt) {
		case 'h' :
			std::cout << helpStream.str();
			return 0;
			break;
		case '!' :
			std::cout << versionStream.str();
			return 0;
			break;
		default:
			// getopt_long has already produced an error message
			return 1;
		}
		opt = getopt_long(argc, argv, OPT_STRING, longOpts, &optionindex);
	}

	// Check that there is still one command line option left: The JSON input file
	if (optionindex >= argc) {
		std::cerr << "Error! No JSON input file given. Use " << argv[0] << " --help "
				  << "for information on how to call SSAGES.\n";
		return 2;
	}

	try{
		// Read in input file
		std::string input = argv[optionindex];
		auto json = JSONLoader().LoadFile(input, world);

		// Build driver.
		auto* driver = Driver::Build(json, world);
		driver->Run();	
	} catch(BuildException& e) {
		if(world.rank() == 0)
			DumpErrorsToConsole(e.GetErrors(), 30);
		MPI_Abort(world, -1);
	} catch(std::exception& e) {
		if(world.rank() == 0)		
			DumpErrorsToConsole({e.what()},	 30);
		MPI_Abort(world, -1);
	} catch(int& k) { 
		std::string err = strerror(k);
		if(world.rank() == 0)		
			DumpErrorsToConsole({"File IO error: " + err}, 30);
		MPI_Abort(world, -1);
	}
	
	return 0;
}

/**
 * This file is part of
 * SSAGES - Suite for Advanced Generalized Ensemble Simulations
 *
 * Copyright 2016 Ben Sikora <bsikora906@gmail.com>
 *                Hythem Sidky <hsidky@nd.edu>
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

#include "../Drivers/Driver.h"
#include "../Grids/Grid.h"
#include "../Constraints/Constraint.h"
#include "../JSON/Serializable.h"
#include "json/json.h"
#include <boost/mpi.hpp>
#include "../JSON/JSONLoader.h"
#include "../JSON/Serializable.h"
#include <exception>
#include <stdexcept>
#include "../Utility/BuildException.h"
#include "../Validator/ArrayRequirement.h"
#include "schema.h"
#include "config.h"

#ifdef ENABLE_GROMACS
#include "Drivers/GromacsDriver.h"
#endif

#ifdef ENABLE_LAMMPS
#include "Drivers/LammpsDriver.h"
#endif

namespace mpi = boost::mpi;
using namespace Json;
namespace SSAGES
{

	//! Simulation class for creating simulation and a driver object.
	/*!
	 * The simulation class pulls together the simulation and the driver to
	 * steer the simulation.
	 *
	 * \ingroup Core
	 */
	class Simulation
	{

	protected:

		//! The MPI world communicator
		boost::mpi::communicator _world;

		//! The local communicator
		boost::mpi::communicator _comm;

		//! A pointer to the driver
		Driver* _MDDriver;

		//! Number of walkers
		int _nwalkers;

		//! MDEngine of choice
		std::string _MDEngine;

		//! Global input file
		std::string _GlobalInput;

		int _ltot; //!< Magic number for message passing
		int _msgw; //!< Magic number for message passing
		int _notw; //!< Magic number for message passing

	public:

		//! Constructor
		/*!
		 * \param world The MPI world communicator to use
		 *
		 * Construct the simulation object.
		 */
		Simulation(boost::mpi::communicator& world) : 
		_world(world), _comm(), _MDDriver(), _nwalkers(1),
		_MDEngine(), _GlobalInput(),
		_ltot(81), _msgw(51), _notw(_ltot - _msgw)
		{

		}

		//! Destructor
		~Simulation()
		{
			delete _MDDriver;
		}

		//! Read the JSON file
		/*!
		 * \param jfile File name of the JSON file
		 * \returns JSON Value containing all input information.
		 *
		 * With this function, the MPI head node will read in the JSON file and
		 * broadcast to other nodes the necessary information. The JSON file
		 * will include the Engine input file name(s).
		 */
		Json::Value ReadJSON(std::string jfile)
		{
			Value root;
			JSONLoader loader;

			if(_world.rank() == 0)
				PrintBoldNotice(" > Validating JSON...", _msgw, _world);

			try{
				root = loader.LoadFile(jfile, _world);
			} catch(std::exception& e) {
				
				DumpErrorsToConsole({e.what()}, _notw);
				_world.abort(-1);
			} catch(int& k) { 
				std::string err = strerror(k);
				DumpErrorsToConsole({"File IO error: " + err}, _notw);
				_world.abort(-1);
			}

				if(_world.rank() == 0)
					std::cout << std::setw(_notw) << std::right << "\033[32mOK!\033[0m\n";

			return root;

		}

		//! Set up the simulation
		/*!
		 * \param json JSON value containing the input information.
		 * \param path Path for JSON path specification.
		 * \return \c True if simulation was set up correctly, else return \c False .
		 *
		 * This function sets up the simulation based on the JSON input
		 * information read in before using Simulation::ReadJSON().
		 */
		bool BuildSimulation(const Json::Value& json, const std::string& path)
		{

			ObjectRequirement validator;
			Value schema;
			Reader reader;

			reader.parse(JsonSchema::Simulation, schema);
			validator.Parse(schema, path);
			
			try
			{
				// Validate inputs.
				validator.Validate(json, path);
				if(validator.HasErrors())
				{
					throw BuildException(validator.GetErrors());
				}
			} catch(BuildException& e) { 	
				DumpErrorsToConsole(e.GetErrors(), _notw);
				return false;
			} catch(std::exception& e) {
				DumpErrorsToConsole({e.what()}, _notw);
				return false;
			}


			_GlobalInput = json.get("inputfile","none").asString();

			return true;
		}

		//! Set up the driver for the simulation
		/*!
		 * \param json JSON value containing input information.
		 * \param path JSON path specification.
		 * \return \c True if Driver was set up correctly, else return \c False .
		 *
		 * Set up the Driver for the Metadynamics simulation.
		 */
		bool BuildDriver(const Json::Value& json, const std::string& path)
		{

			ArrayRequirement validator;
			Value schema;
			Value JsonDrivers;
			Value JsonDriver;
			Reader reader;

			bool success_build = true;

			reader.parse(JsonSchema::driver, schema);
			validator.Parse(schema, path);

			JsonDrivers = json.get("driver",Json::arrayValue);
			try
			{
				// Validate inputs.
				validator.Validate(JsonDrivers, path);
				if(validator.HasErrors())
				{
					throw BuildException(validator.GetErrors());
				}
			} catch(BuildException& e) {
				DumpErrorsToConsole(e.GetErrors(), _notw);
				return false;
			} catch(std::exception& e) {
				DumpErrorsToConsole({e.what()}, _notw);
				return false;
			}

			// Determine walker ID using JSON file and split 
			// communicator accordingly. 
			int i = 0;
			int totalproc = 0;
			int wid = 0;
			for(auto& m : JsonDrivers)
			{
				int currentnumproc = m.get("number processors", 1).asInt();

				if ((int)_world.rank() >= totalproc && (int)_world.rank() < currentnumproc + totalproc)
				{
					wid = i;
					JsonDriver = m;
				}

				totalproc += currentnumproc;
				i++;
			}

			_nwalkers = JsonDrivers.size();
			_comm = _world.split(wid);

			if(_world.size() != totalproc)
			{
				if(_world.rank() == 0)
				{
					std::cerr << "The number of processors for each driver (walker) must sum "
					<< "to the total processors for the mpi call." << std::endl;
					std::cerr<<_world.size()<<" vs "<<totalproc<<std::endl;
				}
				return false;
			}

			// Get the engine. 
			_MDEngine = JsonDriver.get("type", "none").asString();

			// Use input from JSON to determine MDEngine of choice as well as other parameters
			#ifdef ENABLE_LAMMPS
			if(_MDEngine == "LAMMPS")
			{
				LammpsDriver* en = new LammpsDriver(_world, _comm, wid);

				if(!(_MDDriver = static_cast<Driver*>(en)))
				{
						std::cerr << "Unable to dynamic cast engine on node "<<_world.rank()<<" Error occurred" << std::endl;
						success_build = false;		
				}
			}
			else
			#endif
			#ifdef ENABLE_GROMACS
			if(_MDEngine == "Gromacs")
			{
				GromacsDriver* en = new GromacsDriver(_world, _comm, wid);
				if(!(_MDDriver = static_cast<Driver*>(en)))
				{
					std::cout << "Unable to static cast engine on node " << _world.rank() << std::endl;
					success_build = false;
				}
			}
			else
			#endif
			{
				DumpErrorsToConsole({"Unknown MD Engine [" + _MDEngine + "] specified."},_notw);
				success_build = false;
			}

			if (!success_build) {
				// If building the MDEngine failed, don't bother about the rest.
				return success_build;
			}

			// Build the driver
			try{
				_MDDriver->BuildDriver(JsonDriver, path + "/" + std::to_string(wid));
			} catch(BuildException& e) {
		        DumpErrorsToConsole(e.GetErrors(), _notw);
		        success_build = false;
			} catch(std::exception& e) {
		        DumpErrorsToConsole({e.what()}, _notw);
		        success_build = false;
			}

            // Build the CVs
			//if(!json.isMember("CVs") && !JsonDriver.isMember("CVs"))
			//{
			//	DumpErrorsToConsole({"Need global CVs or per driver CVs"},_notw);
			//	success_build = false;
			//}
			//else if(JsonDriver.isMember("CVs"))
			if(JsonDriver.isMember("CVs"))
			{
				if(!BuildCVs(JsonDriver, "#/CVs"))
					success_build = false;
			}
			else
			{
				if(!BuildCVs(json,"#/CVs"))
					success_build = false;
			}

			// Build the method
			//if(!json.isMember("method") && !JsonDriver.isMember("method"))			{
			//	DumpErrorsToConsole({"Need global method or per driver method"},_notw);
			//	success_build = false;
			//}
			//else if (JsonDriver.isMember("method"))
			if (JsonDriver.isMember("method"))
			{
				if(!BuildMethod(JsonDriver, "#/Methods"))
					success_build = false;
			}
			else if (json.isMember("method"))
			{
				if(!BuildMethod(json,"#/Methods"))
					success_build = false;
			}

			// Build the grid if it exists
			if(JsonDriver.isMember("grid"))
			{
				if(!BuildGrid(JsonDriver,"#/Grids"))
					success_build = false;
			}
			else if (json.isMember("grid"))
			{
				if(!BuildGrid(json,"#/Grids"))
					success_build = false;
			}

			if(JsonDriver.isMember("observers"))
			{
				if(!BuildObservers(JsonDriver))
					success_build = false;
			}
			else if (json.isMember("observers"))
			{
				if(!BuildObservers(json))
					success_build = false;
			}

			if(_comm.rank()==0)
			{
				if(success_build)
					std::cout << std::setw(_notw) << std::right << "\033[32mMDEngine " << wid <<  " pass!\033[0m\n";
				else
					std::cout << std::setw(_notw) << std::right << "\033[32mMDEngine " << wid <<  " FAIL!\033[0m\n";
			}

			// Build the constraints if they exist
			if(JsonDriver.isMember("constraints"))
			{
				if(!BuildConstraints(JsonDriver, "#/Constraints"))
					success_build = false;
			}
			else if (json.isMember("constraints"))
			{
				if(!BuildConstraints(json, "#/Constraints"))
					success_build = false;
			}

			if(_comm.rank()==0)
			{
				if(success_build)
					std::cout << std::setw(_notw) << std::right << "\033[32mMDEngine " << wid <<  " pass!\033[0m\n";
				else
					std::cout << std::setw(_notw) << std::right << "\033[32mMDEngine " << wid <<  " FAIL!\033[0m\n";
			}

            return success_build;
		}

		//! Set up CVs
		/*!
		 * \param json JSON value containing input information.
		 * \param path JSON path specification.
		 * \return \c True if CVs have been set up correctly, else return \c False .
		 *
		 * Set up the collective variables (CV) used in the Metadynamics
		 * simulation.
		 */
		bool BuildCVs(const Json::Value& json, const std::string& path)
		{
			// Test if MDDriver has been built
			if (!_MDDriver) {
				throw std::runtime_error("Simulation driver needs to be built before CVs.");
			}

			// Build CV(s).
			try{
				_MDDriver->BuildCVs(json.get("CVs", Json::arrayValue), path);
			} catch(BuildException& e) {
				DumpErrorsToConsole(e.GetErrors(), _notw);
				return false;
			} catch(std::exception& e) {
				DumpErrorsToConsole({e.what()}, _notw);
				return false;
			}

			return true;
		}

		bool BuildConstraints(const Json::Value& json, const std::string& path)
		{
			// Build CV(s).
			try{
				_MDDriver->BuildConstraints(json.get("constraints", Json::arrayValue), path);
			} catch(BuildException& e) {
				DumpErrorsToConsole(e.GetErrors(), _notw);
				return false;
			} catch(std::exception& e) {
				DumpErrorsToConsole({e.what()}, _notw);
				return false;
			}

			return true;
		}

		//! Set up Method
		/*!
		 * \param json JSON value containing input information.
		 * \param path JSON path specification.
		 * \return \c True if Method has been set up correctly, else return \c False .
		 *
		 * Set up the Metadynamics method to be used in the simulation.
		 */
		bool BuildMethod(const Json::Value& json, const std::string& path)
		{
			// Test if MDDriver has been built.
			if (!_MDDriver) {
				throw std::runtime_error("Simulation Driver needs to be built before Methods.");
			}

			// Build method(s).
			try{
				_MDDriver->BuildMethod(json.get("method", Json::objectValue), path);
			} catch(BuildException& e) {
				DumpErrorsToConsole(e.GetErrors(), _notw);
				return false;
			} catch(std::exception& e) {
				DumpErrorsToConsole({e.what()}, _notw);
				return false;
			}

			return true;
		}

		//! Set up observers
		/*!
		 * \param json JSON value containing input information.
		 *
		 * \return \c True if Method has been set up correctly, else return \c False .
		 *
		 * Set up the JSON observer to be used in creating restarts for the simulation.
		 */
		bool BuildObservers(const Json::Value& json)
		{
			// Check that MDDriver has been built.
			if (!_MDDriver) {
				throw std::runtime_error(
					"Simulation Driver needs to be built before Observers."
				);
			}

			// Build method(s).
			try{
				_MDDriver->BuildObservers(json.get("observers", Json::objectValue), _nwalkers);
			} catch(BuildException& e) {
				DumpErrorsToConsole(e.GetErrors(), _notw);
				return false;
			} catch(std::exception& e) {
				DumpErrorsToConsole({e.what()}, _notw);
				return false;
			}

			return true;
		}

		//! Set up grid
		/*!
		 * \param json JSON value containing input information.
		 * \param path JSON path specification.
		 * \return \c True if Grid has been set up sucessfully, else return \c False .
		 *
		 * Set up the grid for the Metadynamics simulation.
		 */
		bool BuildGrid(const Json::Value& json, const std::string& path)
		{
			// Check that MDDriver has been built
			if (!_MDDriver) {
				throw std::runtime_error("Simulation Driver needs to be built before Grid.");
			}

			try{
				_MDDriver->BuildGrid(json, path);
			} catch(BuildException& e) {
				DumpErrorsToConsole(e.GetErrors(), _notw);
				return false;
			} catch(std::exception& e) {
				DumpErrorsToConsole({e.what()}, _notw);
				return false;
			}

			return true;
		}

		//! Read the MD Engine input file
		/*!
		 * Read the input file for the MD engine used in the simulation and
		 * broadcast the necessary information to all MPI instances.
		 */
		void ReadInputFile()
		{
			// Check that MDDriver has been built
			if (!_MDDriver) {
				throw std::runtime_error(
					"Trying to read input file before Simulation Driver was built."
				);
			}

			std::string contents ="";
			std::string localInput = _MDDriver->GetInputFile();

			if(_GlobalInput != "none" && localInput == "none")
			{
				_MDDriver->SetInputFile(_GlobalInput);
				if(_world.rank() == 0)
				{
					std::cout<<"Global input file found, first using: "<<_GlobalInput<<std::endl;
					contents = GetFileContents(_GlobalInput.c_str());
				}
			}

			mpi::broadcast(_world, contents, 0);

			// Each node reads in specified input file
			if(localInput != "none")
			{
				if(_comm.rank() == 0)
				{
					std::cout<<"No/overloaded global input file, node " <<_world.rank()<<" using: "<<localInput<<std::endl;
					contents = GetFileContents(localInput.c_str());
				}
				mpi::broadcast(_comm, contents, 0);
			}

			_MDDriver->ExecuteInputFile(contents);

		}

		//! Set up listeners and hook
		void Finalize()
		{
			if(_MDDriver)
				_MDDriver->Finalize();
		}

		//! Perform simulation run
		void Run()
		{
			// Check that MDDriver has been built
			if(!_MDDriver) {
				throw std::runtime_error(
					"Trying to run simulation before Simulation Driver was built."
				);
			}

			_MDDriver->Run();
		}
	};
}

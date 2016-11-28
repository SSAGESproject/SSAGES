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

#include <stdexcept>
#include <exception>
#include <boost/mpi.hpp>
#include "json/json.h"
#include "schema.h"
#include "config.h"
#include "../Constraints/Constraint.h"
#include "../Drivers/Driver.h"
#include "../Drivers/DriverException.h"
#include "../JSON/Serializable.h"
#include "../JSON/JSONLoader.h"
#include "../JSON/Serializable.h"
#include "../Validator/ArrayRequirement.h"

#ifdef ENABLE_GROMACS
#include "Drivers/GromacsDriver.h"
#endif

#ifdef ENABLE_LAMMPS
#include "Drivers/LammpsDriver.h"
#endif

#ifdef ENABLE_QBOX
#include "Drivers/QboxDriver.h"
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
		boost::mpi::communicator world_;

		//! The local communicator
		boost::mpi::communicator comm_;

		//! A pointer to the driver
		Driver* MDDriver_;

		//! Number of walkers
		int nwalkers_;

		//! MDEngine of choice
		std::string MDEngine_;

		//! Global input file
		std::string GlobalInput_;

		int ltot_; //!< Magic number for message passing
		int msgw_; //!< Magic number for message passing
		int notw_; //!< Magic number for message passing

	public:

		//! Constructor
		/*!
		 * \param world The MPI world communicator to use
		 *
		 * Construct the simulation object.
		 */
		Simulation(boost::mpi::communicator& world) : 
		world_(world), comm_(), MDDriver_(), nwalkers_(1),
		MDEngine_(), GlobalInput_(),
		ltot_(81), msgw_(51), notw_(ltot_ - msgw_)
		{

		}

		//! Destructor
		~Simulation()
		{
			delete MDDriver_;
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

			if(world_.rank() == 0)
				PrintBoldNotice(" > Validating JSON...", msgw_, world_);

			try{
				root = loader.LoadFile(jfile, world_);
			} catch(std::exception& e) {
				
				DumpErrorsToConsole({e.what()}, notw_);
				world_.abort(-1);
			} catch(int& k) { 
				std::string err = strerror(k);
				DumpErrorsToConsole({"File IO error: " + err}, notw_);
				world_.abort(-1);
			}

				if(world_.rank() == 0)
					std::cout << std::setw(notw_) << std::right << "\033[32mOK!\033[0m\n";

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
				DumpErrorsToConsole(e.GetErrors(), notw_);
				return false;
			} catch(std::exception& e) {
				DumpErrorsToConsole({e.what()}, notw_);
				return false;
			}


			GlobalInput_ = json.get("inputfile","none").asString();

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
				DumpErrorsToConsole(e.GetErrors(), notw_);
				return false;
			} catch(std::exception& e) {
				DumpErrorsToConsole({e.what()}, notw_);
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

				if ((int)world_.rank() >= totalproc && (int)world_.rank() < currentnumproc + totalproc)
				{
					wid = i;
					JsonDriver = m;
				}

				totalproc += currentnumproc;
				i++;
			}

			nwalkers_ = JsonDrivers.size();
			comm_ = world_.split(wid);

			if(world_.size() != totalproc)
			{
				if(world_.rank() == 0)
				{
					std::cerr << "The number of processors for each driver (walker) must sum "
					<< "to the total processors for the mpi call." << std::endl;
					std::cerr<<world_.size()<<" vs "<<totalproc<<std::endl;
				}
				return false;
			}

			// Get the engine. 
			MDEngine_ = JsonDriver.get("type", "none").asString();

			// Use input from JSON to determine MDEngine of choice as well as other parameters
			#ifdef ENABLE_LAMMPS
			if(MDEngine_ == "LAMMPS")
			{
				LammpsDriver* en = new LammpsDriver(world_, comm_, wid);

				if(!(MDDriver_ = static_cast<Driver*>(en)))
				{
						std::cerr << "Unable to dynamic cast engine on node "<<world_.rank()<<" Error occurred" << std::endl;
						success_build = false;		
				}
			}
			else
			#endif
			#ifdef ENABLE_GROMACS
			if(MDEngine_ == "Gromacs")
			{
				GromacsDriver* en = new GromacsDriver(world_, comm_, wid);
				if(!(MDDriver_ = static_cast<Driver*>(en)))
				{
					std::cout << "Unable to static cast engine on node " << world_.rank() << std::endl;
					success_build = false;
				}
			}
			else
			#endif
			#ifdef ENABLE_QBOX
			if(MDEngine_ == "Qbox")
			{
				QBoxDriver* en = new QBoxDriver(world_, comm_, wid);
				if(!(MDDriver_ = static_cast<Driver*>(en)))
				{
					std::cout << "Unable to static cast engine on node " << world_.rank() << std::endl;
					success_build = false;
				}
			}
			else
			#endif
			{
				DumpErrorsToConsole({"Unknown MD Engine [" + MDEngine_ + "] specified."},notw_);
				success_build = false;
			}

			if (!success_build) {
				// If building the MDEngine failed, don't bother about the rest.
				return success_build;
			}

			// Build the driver
			try{
				MDDriver_->BuildDriver(JsonDriver, path + "/" + std::to_string(wid));
			} catch(BuildException& e) {
		        DumpErrorsToConsole(e.GetErrors(), notw_);
		        success_build = false;
			} catch(std::exception& e) {
		        DumpErrorsToConsole({e.what()}, notw_);
		        success_build = false;
			}

            // Build the CVs
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

			if(comm_.rank()==0)
			{
				if(success_build)
					std::cout << std::setw(notw_) << std::right << "\033[32mMDEngine " << wid <<  " pass!\033[0m\n";
				else
					std::cout << std::setw(notw_) << std::right << "\033[32mMDEngine " << wid <<  " FAIL!\033[0m\n";
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
			if (!MDDriver_) {
				throw std::runtime_error("Simulation driver needs to be built before CVs.");
			}

			// Build CV(s).
			try{
				MDDriver_->BuildCVs(json.get("CVs", Json::arrayValue), path);
			} catch(BuildException& e) {
				DumpErrorsToConsole(e.GetErrors(), notw_);
				return false;
			} catch(std::exception& e) {
				DumpErrorsToConsole({e.what()}, notw_);
				return false;
			}

			return true;
		}

		//! Set up Constraints.
		/*!
		 * \param json JSON value containing input information.
		 * \param path Path for JSON path specification.
		 * \return \c True if all constraints have been successfully built, else
		 *         return \c False.
		 */
		bool BuildConstraints(const Json::Value& json, const std::string& path)
		{
			// Build CV(s).
			try{
				MDDriver_->BuildConstraints(json.get("constraints", Json::arrayValue), path);
			} catch(BuildException& e) {
				DumpErrorsToConsole(e.GetErrors(), notw_);
				return false;
			} catch(std::exception& e) {
				DumpErrorsToConsole({e.what()}, notw_);
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
			if (!MDDriver_) {
				throw std::runtime_error("Simulation Driver needs to be built before Methods.");
			}

			// Build method(s).
			try{
				MDDriver_->BuildMethod(json.get("method", Json::objectValue), path);
			} catch(BuildException& e) {
				DumpErrorsToConsole(e.GetErrors(), notw_);
				return false;
			} catch(std::exception& e) {
				DumpErrorsToConsole({e.what()}, notw_);
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
			if (!MDDriver_) {
				throw std::runtime_error(
					"Simulation Driver needs to be built before Observers."
				);
			}

			// Build method(s).
			try{
				MDDriver_->BuildObservers(json.get("observers", Json::objectValue), nwalkers_);
			} catch(BuildException& e) {
				DumpErrorsToConsole(e.GetErrors(), notw_);
				return false;
			} catch(std::exception& e) {
				DumpErrorsToConsole({e.what()}, notw_);
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
			if (!MDDriver_) {
				throw std::runtime_error(
					"Trying to read input file before Simulation Driver was built."
				);
			}

			std::string contents ="";
			std::string localInput = MDDriver_->GetInputFile();

			if(GlobalInput_ != "none" && localInput == "none")
			{
				MDDriver_->SetInputFile(GlobalInput_);
				if(world_.rank() == 0)
				{
					std::cout<<"Global input file found, first using: "<<GlobalInput_<<std::endl;
					contents = GetFileContents(GlobalInput_.c_str());
				}
			}

			mpi::broadcast(world_, contents, 0);

			// Each node reads in specified input file
			if(localInput != "none")
			{
				if(comm_.rank() == 0)
				{
					std::cout<<"No/overloaded global input file, node " <<world_.rank()<<" using: "<<localInput<<std::endl;
					contents = GetFileContents(localInput.c_str());
				}
				mpi::broadcast(comm_, contents, 0);
			}

			MDDriver_->ExecuteInputFile(contents);

		}

		//! Set up listeners and hook
		void Finalize()
		{
			if(MDDriver_)
				MDDriver_->Finalize();
		}

		//! Perform simulation run
		void Run()
		{
			// Check that MDDriver has been built
			if(!MDDriver_) {
				throw std::runtime_error(
					"Trying to run simulation before Simulation Driver was built."
				);
			}

			MDDriver_->Run();
		}
	};
}

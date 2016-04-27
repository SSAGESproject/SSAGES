#pragma once

#include "Drivers/Driver.h"
#include "Drivers/LammpsDriver.h"
#include "../JSON/Serializable.h"
#include "json/json.h"
#include <boost/mpi.hpp>
#include "../JSON/JSONLoader.h"
#include "../JSON/Serializable.h"
#include <exception>
#include "../Utility/BuildException.h"
#include "../Validator/ArrayRequirement.h"

namespace mpi = boost::mpi;
using namespace Json;
namespace SSAGES
{
	// Simulation class for creating simulation and a driver object.

	class Simulation : public Serializable
	{

	protected:

		// The MPI world communicator
		boost::mpi::communicator _world;

		// The local communicator 
		boost::mpi::communicator _comm;

		// A pointer to the driver
		Driver* _MDDriver;

		// Number of walkers
		int _nwalkers;

		// MDEngine of choice
		std::string _MDEngine;

		// Global input file
		std::string _GlobalInput;

		// For message parsing
		int _ltot, _msgw, _notw;

	public:

		Simulation(boost::mpi::communicator& world) : 
		_world(world), _MDDriver(), _nwalkers(1),
		_ltot(81), _msgw(51), _notw(_ltot - _msgw)
		{

		}

		~Simulation()
		{
			delete _MDDriver;
		}

		Json::Value ReadJSON(std::string jfile)
		{
			Value root;
			JSONLoader loader;

			// Read in JSON using head node and broadcast to other nodes
			// JSON file will include the Engine input file name(s)

			if(_world.rank() == 0)
			{
				// Parse JSON.
				PrintBoldNotice(" > Validating JSON...", _msgw, _world);

				try{
					root = loader.LoadFile(jfile, _world);
				} catch(std::exception& e) {
					if(_world.rank() == 0)
						DumpErrorsToConsole({e.what()}, _notw);
					_world.abort(-1);
				} catch(int& k) { 
					std::string err = strerror(k);
					
					if(_world.rank() == 0)
						DumpErrorsToConsole({"File IO error: " + err}, _notw);
					_world.abort(-1);
				}

				std::cout << std::setw(_notw) << std::right << "\033[32mOK!\033[0m\n";
			}

			return root;

		}
		void BuildSimulation(const Json::Value& json, const std::string& path)
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
				if(_world.rank() == 0)
					DumpErrorsToConsole(e.GetErrors(), _notw);
				_world.abort(-1);
			} catch(std::exception& e) {
				if(_world.rank() == 0)
					DumpErrorsToConsole({e.what()}, _notw);
				_world.abort(-1);
			}


			_GlobalInput = json.get("inputfile","none").asString();
		}

		void BuildDriver(const Json::Value& json, const std::string& path)
		{

			ArrayRequirement validator;
			Value schema;
			Value JsonDriver;
			Reader reader;

			reader.parse(JsonSchema::Driver, schema);
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
				if(_world.rank() == 0)
					DumpErrorsToConsole(e.GetErrors(), _notw);
				_world.abort(-1);
			} catch(std::exception& e) {
				if(_world.rank() == 0)
					DumpErrorsToConsole({e.what()}, _notw);
				_world.abort(-1);
			}

			// Determine walker ID using JSON file and split 
			// communicator accordingly. 
			int i = 0;
			int totalproc = 0;
			int wid = 0;
			for(auto& m : json)
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

			_nwalkers = json.size();
			_comm = _world.split(wid);

			if(_world.size() != totalproc - 1)
			{
				if(_world.rank() == 0)
					std::cerr << "The number of processors for each driver (walker) must sum "
					<< "to the total processors for the mpi call." << std::endl;
				_world.abort(-1);
			}


			std::string localInput = JsonDriver.get("inputfile", "none").asString(); 
			if(localInput != "none")
				_GlobalInput = localInput;

			// Get the engine. 
			_MDEngine = JsonDriver.get("type", "none").asString();

			// Use input from JSON to determine MDEngine of choice as well as other parameters
			if(_MDEngine == "LAMMPS")
			{
				LammpsDriver* en = new LammpsDriver(_world, _comm, wid);
				_MDDriver = static_cast<Driver*>(en);
			}
			else
			{
				std::cout<<"Errors"<<std::endl;
				throw BuildException({"Unknown MD Engine [" + _MDEngine + "] specified."});
			}

			// Build the driver, cv(s), and method(s)
			_MDDriver->BuildDriver(JsonDriver, path + "/" + std::to_string(wid));
			_MDDriver->BuildCVs();
			_MDDriver->BuildMethod();

			// Read in the global/local input file
			_MDDriver->ReadInputFile();

			// Set up listeners and snapshots
			_MDDriver->Finalize();
		}

		void Run()
		{
			_MDDriver->Run();
		}

		virtual void Serialize(Json::Value& json) const override
		{
			json["number walkers"] = _nwalkers;
			json["MDEngine"] = _MDEngine;
		}
	};
}

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

		// JSON value
		Value _root;

		// Number of walkers
		int _nwalkers;

		// MDEngine of choice
		std::string _MDEngine;

		// For message parsing
		int _ltot, _msgw, _notw;

	public:

		Simulation(boost::mpi::communicator& world) : 
		_world(world), _MDDriver(), _root(), _nwalkers(1),
		_ltot(81), _msgw(51), _notw(_ltot - _msgw)
		{

		}

		~Simulation()
		{
			delete _MDDriver;
		}

		void BuildSimulation(std::string jfile)
		{
			ObjectRequirement validator;
			Reader reader;
			JSONLoader loader;

			// Read in JSON using head node and broadcast to other nodes
			// JSON file will include the Engine input file name(s)

			if(_world.rank() == 0)
			{
				// Parse JSON.
				PrintBoldNotice(" > Validating JSON...", _msgw, _world);

				try{
					_root = loader.LoadFile(jfile, _world);
				} catch(std::exception& e) {
					if(_world.rank() == 0)
						DumpErrorsToConsole({e.what()}, _notw);
					
					_MDDriver = nullptr;
				} catch(int& k) { 
					std::string err = strerror(k);
					
					if(_world.rank() == 0)
						DumpErrorsToConsole({"File IO error: " + err}, _notw);
					
					_MDDriver = nullptr;
				}

				std::cout << std::setw(_notw) << std::right << "\033[32mOK!\033[0m\n";
			}

			// Get requested number of processors per walker and make 
			// sure the number of processors is evenly divisible.
			_nwalkers = _root.get("number procs", 1).asInt();
			if(_world.size() % _nwalkers != 0)
			{
				if(_world.rank() == 0)
					std::cerr << "The number of processors must be evenly "
					<< "divisible by the number of walkers." << std::endl;
				_world.abort(-1);
			}

			std::string path = "#/Drivers";

			// Determine walker ID using integer division and split 
			// communicator accordingly. 
			int wid = (int)_world.rank() / _nwalkers;
			_comm = _world.split(wid);

			reader.parse(JsonSchema::Driver, _root);
			validator.Parse(_root, path);

			// Validate inputs.
			validator.Validate(_root, path);
			if(validator.HasErrors())
				throw BuildException(validator.GetErrors());

			// Get the engine. 
			_MDEngine = _root.get("MDEngine", "none").asString();

			// Use input from JSON to determine MDEngine of choice as well as other parameters
			if(_MDEngine == "LAMMPS")
			{
				LammpsDriver* en = new LammpsDriver(_world, _comm, wid, _root);
				_MDDriver = static_cast<Driver*>(en);
			}
			else
			{
				std::cout<<"Errors"<<std::endl;
				throw BuildException({"Unknown MD Engine [" + _MDEngine + "] specified."});
			}
		}

		void BuildDriver()
		{
			// Build the driver, cv(s), and method(s)
			_MDDriver->BuildDriver();
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

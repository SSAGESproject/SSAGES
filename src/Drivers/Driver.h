#pragma once

#include "../JSON/Serializable.h"
#include "json/json.h"
#include <vector>
#include <boost/mpi.hpp>
#include "../Hook.h"
#include "../Methods/Method.h"
#include "../Snapshot.h"


namespace mpi = boost::mpi;
namespace SSAGES
{
	// Pure abstract driver class for creating driver objects.

	class Driver
	{

	protected:
		boost::mpi::communicator _world, _comm;

		// The node id that this driver belongs to
		const int _wid;

		// The hook that hooks into the MD simulation engine
		Hook* _hook;

		// Driver specific input file (e.g. Lammps.in)
		Json::Value& _json;

		// Input file that driver will use
		std::string _inputfile;

		// The snapshot of your system
		Snapshot* _snapshot;

		//The Method that will be used
		Method* _method;

		// The CVs that will be used
		CVList _CVs;

		// For message parsing
		int _ltot, _msgw, _notw;

	public:

		Driver(boost::mpi::communicator& world, 
			   boost::mpi::communicator& comm,
			   int walkerID,
			   Json::Value& jsonfile,
			   std::string inputfile) : 
		_world(world), _comm(comm), _wid(walkerID), _json(jsonfile), _inputfile(inputfile) {}

		virtual ~Driver()
		{
			for(auto& cv : _CVs)
				delete cv;

			_CVs.clear();

			delete _snapshot;

			delete _method;
		}

		virtual void Run() = 0;

		// Build the driver, which will create the hook and so forth
		virtual void BuildDriver(const Json::Value& json, const std::string path) const = 0;

		// Serialize
		virtual void Serialize(Json::Value& json) const = 0;

		// Read in driver specific input file (e.g. Lammps.in)
		virtual void ExecuteInputFile(std::string contents) const = 0;

		// Build CVs
		void BuildCVs()
		{


		}

		void BuildMethod()
		{

		}

		// Create the snapshot and put all gathered values into the local hook
		void Finalize()
		{
			// Initialize snapshot. 
			_snapshot = new Snapshot(_comm, _wid);

			// Set the hook to snapshot
			_hook->SetSnapshot(_snapshot);

			_hook->AddListener(_method);
			for(auto&cv : _CVs)
				_hook->AddCV(cv);
		}

		virtual boost::mpi::communicator GetWorldComm(){return _world;}

		virtual boost::mpi::communicator GetLocalComm(){return _comm;}

		virtual std::string GetInputFile(){return _inputfile;}

	};
}

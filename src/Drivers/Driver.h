#pragma once

#include "../JSON/Serializable.h"
#include "json/json.h"
#include <vector>
#include <boost/mpi.hpp>
#include "../Hook.h"
#include "../Methods/Method.h"
#include "../Snapshot.h"
#include "../JSON/JSONLoader.h"


namespace mpi = boost::mpi;
using namespace Json;
namespace SSAGES
{
	// Abstract driver class for creating driver objects.

	class Driver : public Serializable
	{

	protected:
		boost::mpi::communicator _world, _comm;

		// The node id that this driver belongs to
		const int _wid;

		// The hook that hooks into the MD simulation engine
		Hook* _hook;

		// The snapshot of your system
		Snapshot* _snapshot;

		//The Method that will be used
		Method* _method{nullptr};

		// The CVs that will be used
		CVList _CVs;

		//local input file
		std::string _inputfile;

		// Random number generators for setting seeds
		std::random_device _rd;
		std::mt19937 _gen;

	public:

		Driver(boost::mpi::communicator& world, 
			   boost::mpi::communicator& comm,
			   int walkerID) : 
		_world(world), _comm(comm), _wid(walkerID),
		_hook(), _snapshot(), _method(), _CVs(),
		_rd(), _gen(_rd())
		 {}

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
		virtual void BuildDriver(const Json::Value& json, const std::string& path) = 0;

		// Serialize
		virtual void Serialize(Json::Value& json) const = 0;

		// Read in driver specific input file (e.g. Lammps.in)
		virtual void ExecuteInputFile(std::string contents) = 0;

		std::string GetInputFile(){return _inputfile;}

		// Build CVs
		void BuildCVs(const Json::Value& json, const std::string& path)
		{
			// // Build CV(s).
			// if(_CVs.size()>0)
			// 	for(auto& cv :_CVs)
			// 		delete cv;
			
			// _CVs.clear();

			CollectiveVariable::BuildCV(json, _CVs, path);
		}

		void BuildMethod(const Json::Value& json, const std::string& path)
		{
			// Build method(s).
			std::cout << static_cast<void*>(nullptr) << std::endl;
			delete _method;
			_method = Method::BuildMethod(json, _world, _comm, path);
		}

		void BuildGrid(const Json::Value& json, const std::string& path)
		{
			// Build the grid.
			_method->BuildGrid(json, path);
		}

		// Create the snapshot and put all gathered values into the local hook
		// Set up listeners as well
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
	};
}

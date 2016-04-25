#pragma once

#include "../JSON/Serializable.h"
#include "json/json.h"
#include <vector>
#include <boost/mpi.hpp>
#include "../Hook.h"
#include "../Methods/Method.h"
#include "../Snapshot.h"
#include "../JSON/JSONLoader.h"
#include <exception>
#include "../Utility/BuildException.h"
#include "FileContents.h"


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
		Method* _method;

		// The CVs that will be used
		CVList _CVs;

		// JSON value code
		Value _root;

		//Global input file
		std::string _inputfile;

		// For message parsing
		int _ltot, _msgw, _notw;

	public:

		Driver(boost::mpi::communicator& world, 
			   boost::mpi::communicator& comm,
			   int walkerID,
			   Value root) : 
		_world(world), _comm(comm), _wid(walkerID),
		_hook(), _snapshot(), _method(), _CVs(), _root(root),
		_ltot(81), _msgw(51), _notw(_ltot - _msgw)
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
		virtual void BuildDriver() = 0;

		// Serialize
		virtual void Serialize(Json::Value& json) const = 0;

		// Read in driver specific input file (e.g. Lammps.in)
		virtual void ExecuteInputFile(std::string contents) = 0;

		boost::mpi::communicator GetWorldComm(){return _world;}

		boost::mpi::communicator GetLocalComm(){return _comm;}

		std::string GetInputFile(){return _inputfile;}

		// Build CVs
		bool BuildCVs()
		{
			// Build CV(s).
			PrintBoldNotice(" > Building CV(s)...", _msgw, _world); 
			try{
				CollectiveVariable::BuildCV(
						_root.get("CVs", Json::arrayValue), 
						_CVs);
			} catch(BuildException& e) {
				DumpErrorsToConsole(e.GetErrors(), _notw);
				return false;
			} catch(std::exception& e) {
				DumpErrorsToConsole({e.what()}, _notw);
				return false;
			}
			return true;
		}

		bool BuildMethod()
		{
			// Build method(s).
			PrintBoldNotice(" > Building method(s)...", _msgw, _world); 
			try{
				_method = Method::BuildMethod(
						_root.get("methods", Json::arrayValue),
						_world, _comm, _wid);
			} catch(BuildException& e) {
				DumpErrorsToConsole(e.GetErrors(), _notw);
				return false;
			} catch(std::exception& e) {
				DumpErrorsToConsole({e.what()}, _notw);
				return false;
			}
			
			return true;
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

		void ReadInputFile()
		{

			_inputfile = _root.get("inputfile", "none").asString();
			
			std::string contents;

			// All nodes get the same input file
			if( _inputfile == "none" && _method->GetInputFile() == "none")
			{
				BuildException e({"No input file defined in global scope or methods scope!"});
				DumpErrorsToConsole(e.GetErrors(), _notw);
			}
			// Each node reads in specified input file
			else if(_method->GetInputFile() != "none")
			{
				if(_comm.rank() == 0)
				{
					std::cout<<"No/overloaded global input file, node " <<_wid<<" using: "<<_method->GetInputFile()<<std::endl;
					contents = GetFileContents(_method->GetInputFile().c_str());
				}
				mpi::broadcast(_comm, contents, 0);
			}
			// All nodes get the same input file
			else
			{
				if(_world.rank() == 0)
					contents = GetFileContents(_inputfile.c_str());
				mpi::broadcast(_world, contents, 0);
			}

			ExecuteInputFile(contents);

		}
	};
}

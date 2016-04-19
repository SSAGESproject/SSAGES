#include "SimBuilder.h"
#include "../JSON/JSONLoader.h"
#include "config.h"
#include <boost/mpi.hpp>


using namespace std;

namespace SSAGES
{
	int DumpErrorsToConsole(const vector<string>& msgs, int notw)
	{
		cout << setw(notw) << right << "\033[1;31mError(s)! See below.\033[0m\n";
		for(auto& msg : msgs)
				cout << " * " << msg << "\n";
		return -1;
	}

	void DumpNoticesToConsole(const vector<string>& msgs, string prefix, int notw)
	{

		boost::mpi::communicator comm;
		if(comm.rank() == 0)
		{
			cout << setw(notw) << right << "\033[32mOK!\033[0m\n";
			if(msgs.size() == 0)
				return;
			
			for(auto& msg : msgs)
				cout << prefix << " * " << msg << "\n";
		}
	}

	void PrintBoldNotice(const string& notice, int msgw)
	{
		if(comm.rank() == 0)
			cout << setw(msgw + 8) << left << "\033[1m" + notice + "\033[0m";
	}

	bool SimBuilder::BuildSimulation(const std::string& filename)
	{
		Json::Value root;
		vector<string> notices;
		JSONLoader loader;

		// Parse JSON.
		PrintBoldNotice(" > Validating JSON...", _msgw);

		try{
			root = loader.LoadFile(filename);
		} catch(std::exception& e) {
			
			if(comm.rank() == 0)
				DumpErrorsToConsole({e.what()}, _notw);
			
			return false;
		} catch(int& k) { 
			std::string err = strerror(k);
			
			if(comm.rank() == 0)
				DumpErrorsToConsole({"File IO error: " + err}, _notw);
			return false;
		}

		if(comm.rank() == 0)
			cout << setw(_notw) << right << "\033[32mOK!\033[0m\n";

		// Build Method.
		PrintBoldNotice(" > Building method...", _msgw); 
		try{
			auto* Method = Method::BuildMethod(root["Method"][0],_world, _comm);
			_worlds.push_back(world);

			// Add world to world manager.
			_wm.AddWorld(world);

			// Print notices.
			notices.push_back("Building world \"" + world->GetStringID() + "\"...");
			auto dim = world->GetHMatrix();
							
		} catch(BuildException& e) {
			DumpErrorsToConsole(e.GetErrors(), _notw);
			return false;
		} catch(exception& e) {
			DumpErrorsToConsole({e.what()}, _notw);
			return false;
		}

		if(_worlds.size() == 0)
		{
			DumpErrorsToConsole({"No worlds have been specified."}, _notw);
			return false;
		}
		DumpNoticesToConsole(notices, "",_notw);
		notices.clear();

		// Build forcefield(s).
		PrintBoldNotice(" > Building forcefield(s)...", _msgw); 
		try{
			ForceField::BuildForceFields(
					root.get("forcefields", Json::arrayValue), 
					&_ffm, 
					_forcefields);
		} catch(BuildException& e) {
			DumpErrorsToConsole(e.GetErrors(), _notw);
			return false;
		} catch(exception& e) {
			DumpErrorsToConsole({e.what()}, _notw);
			return false;
		}

		// Make sure we've created forcefields.
		if(_forcefields.size() == 0)
		{
			DumpErrorsToConsole({"No forcefields have been specified."}, _notw);
			return false;
		}

		notices.push_back("Initialized " +  to_string(_forcefields.size()) + " forcefield(s).");

		DumpNoticesToConsole(notices, "",_notw);
		notices.clear();

		// Build constraints.
		PrintBoldNotice(" > Building constraint(s)...", _msgw);
		try{
			auto ccs = root.get("forcefields", Json::arrayValue).get("constraints", Json::arrayValue);
			Constraint::BuildConstraints(ccs, &_ffm, &_wm, _constraints);
		} catch(BuildException& e) {
			DumpErrorsToConsole(e.GetErrors(), _notw);
			return false;
		} catch(exception& e) {
			DumpErrorsToConsole({e.what()}, _notw);
			return false;
		}

		notices.push_back("Initialized " +  to_string(_constraints.size()) + " constraint(s).");
		DumpNoticesToConsole(notices, "",_notw);
		notices.clear();

		// Build move(s).
		PrintBoldNotice(" > Building move(s)...", _msgw); 
		try{
			Move::BuildMoves(root.get("moves", Json::arrayValue), &_mm, &_wm, _moves);
		} catch(BuildException& e) {
			DumpErrorsToConsole(e.GetErrors(), _notw);
			return false;
		} catch(exception& e) {
			DumpErrorsToConsole({e.what()}, _notw);
			return false;
		}

		// Make sure we've created moves.
		if(_moves.size() == 0)
		{
			DumpErrorsToConsole({"No moves have been specified."}, _notw);
			return false;
		}

		for(auto& m : _moves)
			notices.push_back("Initialized " + m->GetName() + " move.");

		DumpNoticesToConsole(notices, "",_notw);
		notices.clear();

		// Only build observers on master rank.
		#ifdef MULTI_WALKER
		if(comm.rank() == 0)
		{
		#endif

		// Build observers.
		PrintBoldNotice(" > Building observer(s)...", _msgw); 
		try{
			SimObserver::BuildObservers(root.get("observers", Json::arrayValue), _observers);
		} catch(BuildException& e) {
			DumpErrorsToConsole(e.GetErrors(), _notw);
			return false;
		} catch(exception& e) {
			DumpErrorsToConsole({e.what()}, _notw);
			return false;
		}

		// Make sure we've created at least one observer.
		if(_observers.size() == 0)
		{
			DumpErrorsToConsole({"No observers have been specified."}, _notw);
			return false;
		}

		for(auto& o : _observers)
		{
			notices.push_back("Initialized " + o->GetName() + " observer.");
			notices.push_back("Set sampling frequency to " + 
				to_string(o->GetFrequency()) +  " sweeps.");
		}

		DumpNoticesToConsole(notices, "",_notw);
		notices.clear();

		// End build walkers on master rank.
		#ifdef MULTI_WALKER
		}
		#endif

		if(root.isMember("histogram")) 
		{
			PrintBoldNotice(" > Building histogram...", _msgw); 
			try{
				_hist = Histogram::BuildHistogram(root["histogram"]);
			} catch(BuildException& e) {
				DumpErrorsToConsole(e.GetErrors(), _notw);
				return false;
			} catch(exception& e) {
				DumpErrorsToConsole({e.what()}, _notw);
				return false;
			}

			notices.push_back("Initialized histogram range [" + 
				to_string(_hist->GetMinimum()) + ", " + 
				to_string(_hist->GetMaximum()) + "].");
			notices.push_back("Created " + to_string(_hist->GetBinCount()) + " bins.");

			DumpNoticesToConsole(notices, "",_notw);
			notices.clear();
		}

		if(root.isMember("orderparameter"))
		{
			PrintBoldNotice(" > Building order parameter...", _msgw);
			try{
				_orderp = DOSOrderParameter::Build(root["orderparameter"], _hist, &_wm);
			} catch(BuildException& e) {
				DumpErrorsToConsole(e.GetErrors(), _notw);
				return false;
			} catch(exception& e) {
				DumpErrorsToConsole({e.what()}, _notw);
				return false;
			}

			DumpNoticesToConsole(notices, "",_notw);
			notices.clear();
		}

		// Build simulation.
		PrintBoldNotice(" > Building simulation...", _msgw);
		try{
			_sim = Simulation::BuildSimulation(root, &_wm, &_ffm, &_mm, _orderp, _hist);

			// Add observers.
			for(auto& o  : _observers)
				_sim->AddObserver(o);

		} catch(BuildException& e) {
			DumpErrorsToConsole(e.GetErrors(), _notw);
			return false;
		} catch(exception& e) {
			DumpErrorsToConsole({e.what()}, _notw);
			return false;
		}

		DumpNoticesToConsole(notices, "",_notw);
		notices.clear();

		return true;
	}
}
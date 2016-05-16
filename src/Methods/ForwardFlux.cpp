#include "ForwardFlux.h"
#include <iostream>
#include "../FileContents.h"
#include <random>

// This method involves a lot of bookkeeping. Typically the world node
// will hold gather all needed information and pass it along as it occurs.

namespace SSAGES
{
	void ForwardFlux::PreSimulation(Snapshot* snap, const CVList& cvs)
	{
		_indexfile.open(_indexfilename.c_str());
		_resultsfile.open(_resultsfilename.c_str());

		_currentnode = AtInterface(cvs);
	}

	void ForwardFlux::PostIntegration(Snapshot* snapshot, const CVList& cvs)
	{
		switch(_restart)
		{
			case NEW:
			{
				SetUpNewLibrary(snapshot, cvs);
				return;
			}
			case LIBRARY:
			{
				_currenthash = 10000*_world.rank();
				if(_world.rank() == 0 && _currentstartingpoint != 0)
				{
					for(auto& value : _successes)
						_resultsfile<<value<<" ";
					_resultsfile<<"\n";
				}

				_indexcontents = "";

				for(auto& s : _localsuccesses)
					s = 0;

				for(auto& s : _successes)
					s = 0;

				_currentnode = 0;

				std::vector<std::vector<std::string> > tmp;
				if(!(ExtractInterfaceIndices(0, _librarycontents, tmp)))
				{
					std::cout<< "Could not locate any files";
					std::cout<< " at first interface!"<<std::endl;
					std::cout<< _librarycontents<<std::endl;
					PostSimulation(snapshot, cvs);
					_world.abort(0);
				}

				if(_currentstartingpoint >= tmp.size())
				{
					if(_world.rank()==0)
						std::cout<<"Ending Simulation..."<<std::endl;
					PostSimulation(snapshot, cvs);
					_world.abort(0);
				}

				std::string dumpfilecontents;
				_shootingconfigfile = tmp[_currentstartingpoint][1];
				if(_world.rank() == 0)
					dumpfilecontents = GetFileContents(_shootingconfigfile.c_str());
				mpi::broadcast(_world, dumpfilecontents, 0);

				ReadConfiguration(snapshot, dumpfilecontents);
				_restart = NONE;
				_currentstartingpoint++;
				return;
			}
			case RESTART:
			{
				if(_world.rank()==0)
					std::cout<<"Running from restart "<<_restartfilename<<std::endl;
				std::string dumpfilecontents;
				_shootingconfigfile = _restartfilename;
				if(_world.rank() == 0)
					dumpfilecontents = GetFileContents(_shootingconfigfile.c_str());
				mpi::broadcast(_world, dumpfilecontents, 0);

				ReadConfiguration(snapshot, dumpfilecontents);
				_restart = NONE;
				return;
			}

			case NEWCONFIG:
			{
				_currenthash = 10000*_world.rank();
				mpi::all_reduce(_world, _indexcontents, _globalcontents, std::plus<std::string>());

				std::string dumpfilecontents;
				_shootingconfigfile = PickConfiguration(_currentnode, _globalcontents);

				if(_world.rank() == 0)
					dumpfilecontents = GetFileContents(_shootingconfigfile.c_str());
				
				mpi::broadcast(_world, dumpfilecontents, 0);
				ReadConfiguration(snapshot, dumpfilecontents);
				_restart = NONE;

				return;
			}
			default:
			{
				break;
			}
		}

		// Locate the interface you are at and check if:
		// Returned to origin or at next interface
		int atinter = AtInterface(cvs);
		if(atinter == _currentnode + 1 || atinter == 0)
		{
			_currentnode++;			
			// Check if you made it to the next one!
			if(atinter == _currentnode)
			{
				auto& ID = snapshot->GetSnapshotID();
				ID = "dump_"+std::to_string(_currentnode)+"_"+std::to_string(_currenthash)+".dump";
				_currenthash++;
				if(_comm.rank()==0)
					WriteConfiguration(snapshot);
				_localsuccesses[_currentnode]++;
			}

			mpi::all_reduce(_world, _localsuccesses[_currentnode], _successes[_currentnode], std::plus<int>());

			if(_successes[_currentnode] > 0)
			{
				_restart = NEWCONFIG;

				// If you have reached the final interface you are "done" with that path
				if(_currentnode >= _centers.size()-1)
				{
					std::cout<< "Found finishing configuration! "<<std::endl;
					_restart = LIBRARY;
				}
			}
			else
			{
				_restart = LIBRARY;
			}
		}
	}

	void ForwardFlux::PostSimulation(Snapshot*, const CVList&)
	{
		// Calculate probabilites for each interface.
		// Write probabilities and paths to file
		mpi::all_reduce(_world, mpi::inplace(_globalcontents), std::plus<std::string>());
		//Close local and global files
		if(_world.rank() == 0)
		{
			_indexfile<<_indexcontents<<std::endl;
			_resultsfile<<"flux in: "<<_fluxin<<std::endl;
			_resultsfile<<"flux out: "<<_fluxout<<std::endl;

			_indexfile.close();
			_resultsfile.close();
		}
	}

	// Setting up new run, so setup new starting configurations at the first interface
	void ForwardFlux::SetUpNewLibrary(Snapshot* snapshot, const CVList& cvs)
	{
		// Get CV values check if at next interface, if so store configuration
		int interface = AtInterface(cvs);
		_shootingconfigfile = "Origin";

		// Flux out of A
		if(interface == 1 && _currentnode == 0)
		{
			auto& ID = snapshot->GetSnapshotID();
			ID = "dump_"+std::to_string(interface)+"_"+std::to_string(_currenthash)+".dump";
			_currenthash++;
			if(_comm.rank()==0)
				WriteConfiguration(snapshot);

			_localsuccesses[_currentnode]++;
			_fluxout++;
		}
		// Flux back in towards A
		else if(interface == 0 && _currentnode == 1)
			_fluxin++;

		_currentnode = interface;

		std::vector<std::vector<std::string> > TempLibrary;
		ExtractInterfaceIndices(0, _indexcontents, TempLibrary);

		if(TempLibrary.size() >= _requiredconfigs && _currentnode == 0)
		{
			_restart = LIBRARY;
			mpi::all_reduce(_world, _indexcontents, _librarycontents, std::plus<std::string>());
			_currentstartingpoint = 0;
		}
	}

	void ForwardFlux::WriteConfiguration(Snapshot* snapshot)
	{
		const auto& positions = snapshot->GetPositions();
		const auto& velocities = snapshot->GetVelocities();
		const auto& atomID = snapshot->GetAtomIDs();
		const auto& dumpfilename = snapshot->GetSnapshotID();

		// Write the dump file out
		std::ofstream dumpfile;
 		dumpfile.open(dumpfilename.c_str());

 		for(size_t i = 0; i< atomID.size(); i++)
 		{
 			dumpfile<<atomID[i]<<" ";
 			dumpfile<<positions[i][0]<<" "<<positions[i][1]<<" "<<positions[i][2]<<" ";
 			dumpfile<<velocities[i][0]<<" "<<velocities[i][1]<<" "<<velocities[i][2]<<std::endl;
		}

		std::vector<std::string> tmpstr;
		tmpstr.push_back(std::to_string(_currentnode));
		tmpstr.push_back(dumpfilename);
 		tmpstr.push_back(_shootingconfigfile);

 		// Update index file of new configuration
 		_indexcontents += tmpstr[0]+" "+tmpstr[1]+" "+tmpstr[2]+"\n";
 		dumpfile.close();
	}

	// Extract all indices for a given interface and contetns. 
	// Return false if couldnt locate anything at a given interface
	bool ForwardFlux::ExtractInterfaceIndices(int interface, const std::string& contents,
											 std::vector<std::vector<std::string> >& InterfaceIndices)
	{
		//Extract configuration indices for a given interface int
		std::istringstream f(contents);
		std::string line;
		while (std::getline(f, line))
		{
			std::string buf; // Have a buffer string
			std::stringstream ss(line); // Insert the string into a stream
			std::vector<std::string> tokens; // Create vector to hold our words

			while (ss >> buf)
			    tokens.push_back(buf);

			if(std::stoi(tokens[0]) == interface)
				InterfaceIndices.push_back(tokens);
		}

		if(InterfaceIndices.size() == 0)
			return false;

		return true;
	}

	void ForwardFlux::ReadConfiguration(Snapshot* snapshot, const std::string& FileContents)
	{
		auto& positions = snapshot->GetPositions();
		auto& velocities = snapshot->GetVelocities();
		auto& atomID = snapshot->GetAtomIDs();
		auto& forces = snapshot->GetForces();
		auto& ID = snapshot->GetSnapshotID();

		ID = _shootingconfigfile;

		//Extract _currentconfig information
		std::istringstream f(FileContents);
		std::string line;
		while (std::getline(f, line))
		{
			int atomindex = -1;
			std::string buf; // Have a buffer string
			std::stringstream ss(line); // Insert the string into a stream
			std::vector<std::string> tokens; // Create vector to hold our words

			while (ss >> buf)
			    tokens.push_back(buf);

			if(tokens.size() != 7)
			{
				std::cout<<"error, incorrect line format in "<<_shootingconfigfile<<" on line: "<<std::endl;
				std::cout<<line<<std::endl;
				_world.abort(-1);	
			}

			for(size_t i=0; i < atomID.size(); i++)
			{
				if(atomID[i] == std::stoi(tokens[0]))
					atomindex = i;
			}

			if(atomindex < 0)
			{
				std::cout<<"error, could not locate atomID "<<tokens[0]<<" from dumpfile"<<std::endl;
				_world.abort(-1);
			}

			positions[atomindex][0] = std::stod(tokens[1]);
			positions[atomindex][1] = std::stod(tokens[2]);
			positions[atomindex][2] = std::stod(tokens[3]);
			velocities[atomindex][0] = std::stod(tokens[4]);
			velocities[atomindex][1] = std::stod(tokens[5]);
			velocities[atomindex][2] = std::stod(tokens[6]);

			for(auto& force : forces)
				for(auto& xyz : force)
					xyz = 0;
		}
	}

	// Pick a random configuration
	std::string ForwardFlux::PickConfiguration(int interface, const std::string& contents)
	{
		std::string configfilename;
		if(_world.rank() == 0)
		{
			std::vector<std::vector<std::string> > files; 
			if(!(ExtractInterfaceIndices(interface, contents, files)))
			{
				std::cout<< "Could not locate any files at interface ";
				std::cout<< interface<<" in PickCOnfiguration!"<<std::endl;
				std::cout<< contents<<std::endl;
				_world.abort(-1);
			}

			std::uniform_int_distribution<> dis(0, files.size()-1);
			int configfile = dis(_gen);
			configfilename = files[configfile][1];
		}

		mpi::broadcast(_world, configfilename, 0);

		return configfilename;
	}
}

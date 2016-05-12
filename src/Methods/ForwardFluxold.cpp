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
		// Before the simulation begins:
		// Open up files and read if needed
		if(_world.rank() == 0)
		{
			switch(_restart)
			{
				case LIBRARY:
				{
					_globalcontents = GetFileContents(_indexfilename.c_str());
					_resultscontents = GetFileContents(_resultsfilename.c_str());

					_indexfile.open(_indexfilename.c_str(), std::ofstream::out | std::ofstream::app);
					_resultsfile.open(_resultsfilename.c_str(), std::ofstream::out | std::ofstream::app);
					break;
				}
				case NEW:
				{
					_indexfile.open(_indexfilename.c_str());
					_resultsfile.open(_resultsfilename.c_str());
					break;
				}
				default:
				{
					std::cout<<"Must be LIBRARY or NEW!"<<std::endl;
					_world.abort(0);
					break;
				}
			}
		}

		switch(_restart)
		{
			// Restarting from library
			case LIBRARY:
			{
				if(!SetUpRestartRun(snap))
				{
					if(_world.rank() == 0)
					{
						_indexfile.close();
						_resultsfile.close();
						_indexfile.open(_indexfilename.c_str());
						_resultsfile.open(_indexfilename.c_str());
					}
					_restart = NEW;
					CleanUp();
				}
				break;
			}
			// New simulation
			case NEW:
			{
				CleanUp();
				break;
			}
			default:
			{
				std::cout<<"Error starting forward flux"<<std::endl;
				_world.abort(0);
			}
		}
	}

	void ForwardFlux::PostIntegration(Snapshot* snapshot, const CVList& cvs)
	{
		std::vector<std::vector<std::string> > tmp;
		switch(_restart)
		{
			// Create a new library of starting configs
			case NEW:
			{
				SetUpNewRun(snapshot, cvs);
				return;
			}
			// Run from starting configs
			case LIBRARY:
			{
				_currentstartingpoint++;
				if(!(ExtractInterfaceIndices(1, _librarycontents, tmp)))
				{
					std::cout<< "Could not locate any files";
					std::cout<< " at first interface!"<<std::endl;
					std::cout<< _librarycontents<<std::endl;
					_world.abort(0);
				}

				if(_currentstartingpoint >= tmp.size())
				{
					if(_world.rank()==0)
						std::cout<<"Ending Simulation..."<<std::endl;
					PostSimulation(snapshot, cvs);
					_world.abort(0);
				}

				//ReadConfiguration(snapshot, &_libraryconfigs[_currentstartingpoint]);
				ReadConfiguration(snapshot, tmp[_currentstartingpoint][1]);

				_restart = NONE;
				_indexcontents = "";
				_currentshot = 0;
				_currentnode = 1;
				break;
			}
			// Run from a new configuration on an interface
			case NEWCONFIG:
			{
				PickConfiguration();
				ReadConfiguration(snapshot);
				_indexcontents = "";
				_currentshot = 0;
				_restart = NONE;
				break;
			}
			// Shoot another trajectory using same walker from old configuration
			case OLDCONFIG:
			{
				ReadConfiguration(snapshot, _currentconfig);
				_restart = NONE;
				break;
			}
			default:
			{
				break;
			}
		}

		int atinter = AtInterface(cvs);

		// Check if at interface or back to start
		if(atinter == _currentnode + 1 || atinter == 0)
		{
			_currentshot++;
			_currentnode++;

			// Check if you made it to the next one!
			if(atinter == _currentnode)
			{
				_localsuccesses[_currentnode]++;
				StoreConfiguration(snapshot);
			}

			// Make sure you made all the shots you need to make for this walker
			if(_currentshot < _numshots)
			{
				_currentnode--;
				_restart = OLDCONFIG;
				return;
			}

			mpi::all_reduce(_world, _localsuccesses[_currentnode], _successes[_currentnode], std::plus<int>());
			
			// Did at least one walker make it?
			if(_successes[_currentnode] > 0)
			{
				mpi::all_reduce(_world, mpi::inplace(_indexcontents), std::plus<std::string>());
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
				// Didn't make it so write out this path
				for(size_t k = 0; k<_successes.size();k++)
					_paths[_currentstartingpoint][k] = _successes[k];

				for(auto& l : _localsuccesses)
					l = 0;

				if(_world.rank()==0)
				{
					std::cout<<"Did not advance past "<<_currentnode;
					std::cout<<". Attempting from config "<<_currentstartingpoint<<" on interface 1."<<std::endl;
					std::cout<<std::flush;
					std::cout<<"Node location ";
					for(auto& center : _centers[_currentnode])
						std::cout<<center<<" ";
					std::cout<<std::endl;
				}

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
			_indexfile << _globalcontents<<std::endl;
			_resultsfile<<"flux in: "<<_fluxin<<std::endl;
			_resultsfile<<"flux out: "<<_fluxout<<std::endl;
			_resultsfile<<"Shots: "<<_numshots<<std::endl;
			for(auto& p : _paths)
			{
				for(size_t i = 0; i < p.size(); i++)
					_resultsfile<<p[i]<<" ";
				_resultsfile<<std::endl;
			}

			_indexfile.close();
			_resultsfile.close();
		}
	}

	// Extract all indices for a given interface and contetns. 
	// Return false if couldnt locate anything at a given interface
	bool ForwardFlux::ExtractInterfaceIndices(int interface, const std::string contents,
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

	void ForwardFlux::ReadConfiguration(Snapshot* snapshot, std::string dumpfilename)
	{
		std::string dumpfilecontents;

		_shootingconfigfile = dumpfilename;

		if(_world.rank() == 0)
			dumpfilecontents = GetFileContents(dumpfilename.c_str());

		mpi::broadcast(_world, dumpfilecontents, 0);

		auto& positions = snapshot->GetPositions();
		auto& velocities = snapshot->GetVelocities();
		auto& atomID = snapshot->GetAtomIDs();
		auto& forces = snapshot->GetForces();

		//Extract dump file information
		std::istringstream f(dumpfilecontents);
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
				std::cout<<"error, incorrect line format in "<<dumpfilename<<" on line: "<<std::endl;
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

	void ForwardFlux::ReadConfiguration(Snapshot* snapshot)
	{
		auto& positions = snapshot->GetPositions();
		auto& velocities = snapshot->GetVelocities();
		auto& atomID = snapshot->GetAtomIDs();
		auto& forces = snapshot->GetForces();
		auto& ID = snapshot->GetSnapshotID();

		//Extract _currentconfig information
		std::istringstream f(_currentconfig);
		std::string name;
		std::getline(f,name);

		_shootingconfigfile = name;
		ID = name;
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
				std::cout<<"error, incorrect line format in "<<name<<" on line: "<<std::endl;
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
	std::string ForwardFlux::PickConfiguration(int interface, std::string contents)
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

	// Pick a random configuration
	void ForwardFlux::PickConfiguration()
	{
		_currentconfig = "";
		int NumConfig = _dumpconfigs.size();
		std::vector<int> Configs;
		// Gather total number of configs
		mpi::all_gather(_world, NumConfig, Configs);
		mpi::reduce(_world, NumConfig, std::plus<int>(), 0);

		// Pick a random config
		if(_world.rank() == 0)
		{
			std::uniform_int_distribution<> dis(0, NumConfig-1);
			NumConfig = dis(_gen);
		}

		mpi::broadcast(_world, NumConfig, 0);

		int currentrank = -1;
		int sum = 0;

		while(sum < NumConfig)
		{
			currentrank++;
			sum += Configs[currentrank];
		}

		// Serrialize config
		if(_world.rank() == currentrank)
		{
			const auto& positions = _dumpconfigs[sum - NumConfig].GetPositions();
			const auto& velocities = _dumpconfigs[sum - NumConfig].GetVelocities();
			const auto& atomID = _dumpconfigs[sum - NumConfig].GetAtomIDs();
			const auto& snapshotID = _dumpconfigs[sum - NumConfig].GetSnapshotID();

			std::ostringstream oss;
			// Write the dump file out
			oss<<snapshotID<<std::endl;
	 		for(size_t i = 0; i< atomID.size(); i++)
	 		{
	 			oss<<atomID[i]<<" ";
	 			oss<<positions[i][0]<<" "<<positions[i][1]<<" "<<positions[i][2]<<" ";
	 			oss<<velocities[i][0]<<" "<<velocities[i][1]<<" "<<velocities[i][2]<<std::endl;
			}

			_currentconfig = oss.str();
		}

		mpi::broadcast(_world, _currentconfig, currentrank);
	}

	// Setting up new run, so setup new starting configurations at the first interface
	void ForwardFlux::SetUpNewRun(Snapshot* snapshot, const CVList& cvs)
	{
		// Get CV values check if at next interface, if so store configuration
		int interface = AtInterface(cvs);
		std::vector<std::vector<std::string> > TempLibrary;

		// Flux out of A
		if(interface == 1 && _currentnode == 0)
		{
			_currentnode = interface;
			StoreConfiguration(snapshot);
			WriteConfiguration(snapshot);
			_localsuccesses[_currentnode]++;
			_fluxout++;
		}
		// Flux back in towards A
		else if(interface == 0 && _currentnode == 1)
		{
			_currentnode = interface;
			_fluxin++;
		}

		// Get all starting configs for this walker
		ExtractInterfaceIndices(1, _indexcontents, TempLibrary);

		// See if found the required starting configs
		if(TempLibrary.size() >= _requiredconfigs && _currentnode == 0)
		{
			mpi::all_reduce(_world, _indexcontents, _librarycontents, std::plus<std::string>());
			
			_restart = LIBRARY;
			_indexcontents="";

			int temp=0;
			if(_comm.rank() == 0)
				temp = TempLibrary.size();

			mpi::all_reduce(_world,mpi::inplace(temp),std::plus<int>());
			_paths.resize(temp);
			for(auto& path : _paths)
				path.resize(_successes.size());

			//Start at first location
			_currentstartingpoint = -1;
			_currentnode = 1;
		}
	}

	// Setup for a restart run, currently not implemented
	bool ForwardFlux::SetUpRestartRun(Snapshot* snapshot)
	{
		std::vector<std::vector<std::string> > tmp;
		std::vector<std::vector<std::string> > tmps;
		if(!(ExtractInterfaceIndices(_currentnode, _globalcontents, tmp)))
		{
			std::cout<< "Could not locate any files";
			std::cout<< " reading from library, starting new run!"<<std::endl;
			std::cout<<"Library"<< _globalcontents<<std::endl;
			std::cout<<"interface" << _currentnode<<std::endl;
			return false;
		}

		ReadConfiguration(snapshot, PickConfiguration(_currentnode, _globalcontents));
		_restart = NONE;

		for (auto& success : _successes)
			success = 0;

		_indexcontents = "";
		_librarycontents = "";
		ExtractInterfaceIndices(1, _globalcontents, tmps);
		for(auto& t : tmps)
			_librarycontents += t[0]+" "+t[1]+" "+t[2]+"\n";
		_resultscontents = "";

		_currenthash = 100000*_world.rank();
		_currentstartingpoint = -1;

		_fluxout = _fluxin = 0;
		return true;
	}

	void ForwardFlux::CleanUp()
	{
		for (auto& success : _successes)
			success = 0;

		_indexcontents = "";
		_globalcontents = "";
		_librarycontents = "";
		_resultscontents = "";
		_dumpconfigs.clear();
		_currentconfig.clear();

		_currentnode = 0;
		_currenthash = 100000*_world.rank();
		_currentstartingpoint = -1;

		_fluxout = _fluxin = 0;

	}

	void ForwardFlux::WriteConfiguration(Snapshot* snapshot)
	{

		std::string dumpfilename = "dump_"+std::to_string(_currentnode)+"_"+std::to_string(_currenthash)+".dump";

		if(_comm.rank() == 0)
		{
			const auto& positions = snapshot->GetPositions();
			const auto& velocities = snapshot->GetVelocities();
			const auto& atomID = snapshot->GetAtomIDs();

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

	 		// Update starting library 
	 		if(_currentnode  == 1)
	 			tmpstr.push_back("Origin");		 			
	 		else
	 			tmpstr.push_back(_shootingconfigfile);

	 		// Update index file of new configuration
	 		_indexcontents += tmpstr[0]+" "+tmpstr[1]+" "+tmpstr[2]+"\n";
	 		_globalcontents += tmpstr[0]+" "+tmpstr[1]+" "+tmpstr[2]+"\n";
	 		dumpfile.close();
		}
	}
}
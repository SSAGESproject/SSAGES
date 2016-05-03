#include "ForwardFlux.h"
#include <iostream>
#include "../FileContents.h"
#include <random>

namespace SSAGES
{
	void ForwardFlux::PreSimulation(Snapshot*, const CVList& cvs)
	{	
		// Before the simulation begins:
		// Open up global file for indexing successful configurations
		if(_world.rank() == 0)
		{
			if(_newrun)
			{
		 		_indexfile.open(_indexfilename.c_str());
	 		 	_resultsfile.open(_resultsfilename.c_str());
			}
			else
			{
				_indexcontents = GetFileContents(_indexfilename.c_str());
				_resultscontents = GetFileContents(_resultsfilename.c_str());

				_indexfile.open(_indexfilename.c_str(), std::ofstream::out | std::ofstream::app);
				_resultsfile.open(_resultsfilename.c_str(), std::ofstream::out | std::ofstream::app);
		 	}
		}

		if(!_newrun)
		{
			mpi::broadcast(_world, _indexcontents, 0);
			mpi::broadcast(_world, _resultscontents, 0);
		}

		if(!_newrun)
		{
			if(ExtractInterfaceIndices(_currentinterface, _startinglibrary))
			{
				if(_world.rank() == 0)
				{
					std::cout << "Could not locate an interfaces at " << _currentinterface << " in ";
					std::cout << _indexfilename << "! New starting configurations will be generated." << std::endl;
				}
				
				// Clear and set private variables as if a new run is occuring
				ClearFiles();
			}
		}
	}

	void ForwardFlux::PostIntegration(Snapshot* snapshot, const CVList& cvs)
	{
		if(_newrun)
		{
			SetUpNewRun(snapshot, cvs);
			return;
		}

		if(_restartfromlibrary)
		{
			ReadConfiguration(snapshot, _startinglibrary[_currentstartingpoint][1]);
			_restartfromlibrary = false;
			return;
		}
		else if(_restartfrominterface)
		{
			ReadConfiguration(snapshot, PickConfiguration(_currentinterface));
			_restartfrominterface = false;
			return;
		}

		int atinter = AtInterface(cvs);

		if(atinter == _currentinterface + 1 || atinter == 0)
		{
			if(atinter == _currentinterface + 1)
			{
				_currentinterface++;
				_successes[_currentinterface]++;
				WriteConfiguration(snapshot);
			}

			MPI_Barrier(_world);
			mpi::all_reduce(_world, mpi::inplace(&_successes.front()), _successes.size(), std::plus<int>());

			_restartfromlibrary = true;
			for(auto& success : _successes)
				if(success > 0)
				{
					_restartfromlibrary = false;
					_restartfrominterface = true;
				}

			// If you need to start over, start from next starting library config
			if(_restartfromlibrary == true)
				_currentstartingpoint++;

			// If you have exhausted initial starting configs, start a new run
			if(_currentstartingpoint >= _startinglibrary.size())
			{
				std::cout<<"Could not locate transition from A->B with generated starting library!"<<std::endl;
				_world.abort(0);
			}

			if(_currentinterface >= _centers.size()-1)
			{
				std::cout<< "Found finishing configuration! "<<std::endl;
				PostSimulation(snapshot, cvs);
				_world.abort(0);
			}
		}
	}

	void ForwardFlux::PostSimulation(Snapshot*, const CVList&)
	{
		// Calculate probabilites for each interface.
		// Write probabilities and paths to file

		//Close local and global files
		if(_world.rank() == 0)
		{
			_indexfile.close();
			_resultsfile.close();
		}
	}

	// Extract all indices for a given interface. 
	// Return true if couldnt locate anything at a given interface
	bool ForwardFlux::ExtractInterfaceIndices(int interface, std::vector<std::vector<std::string> >& InterfaceIndices)
	{
		//Extract configuration indices for a given interface int
		std::istringstream f(_indexcontents);
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
			return true;

		return false;
	}

	int ForwardFlux::AtInterface(const CVList& cvs)
	{
		std::vector<double> dists;
		dists.resize(_centers[0].size());

		// Record the difference between all cvs and all nodes
		for (size_t i = 0; i < _centers.size(); i++)
		{
			dists[i] = 0;
			for(size_t j = 0; j < cvs.size(); j++)
				dists[i]+=(cvs[j]->GetValue() - _centers[i][j])*(cvs[j]->GetValue() - _centers[i][j]);
		}

		return (std::min_element(dists.begin(), dists.end()) - dists.begin());
	}

	void ForwardFlux::WriteConfiguration(Snapshot* snapshot)
	{
		const auto& positions = snapshot->GetPositions();
		const auto& velocities = snapshot->GetVelocities();
		const auto& atomID = snapshot->GetAtomIDs();

		if(_comm.rank() == 0)
		{	
			// Write the dump file out
			std::string dumpfilename = "dump_"+std::to_string(_currentinterface)+"_"+std::to_string(_currenthash)+".dump";
			std::ofstream dumpfile;
	 		dumpfile.open(dumpfilename.c_str());

	 		for(size_t i = 0; i< atomID.size(); i++)
	 		{
	 			dumpfile<<atomID[i]<<" ";
	 			dumpfile<<positions[i][0]<<" "<<positions[i][1]<<" "<<positions[i][2]<<" ";
	 			dumpfile<<velocities[i][0]<<" "<<velocities[i][1]<<" "<<velocities[i][2]<<std::endl;
			}

 			std::vector<std::string> tmpstr;
 			tmpstr.push_back(std::to_string(_currentinterface));
 			tmpstr.push_back(dumpfilename);
	 		// Update starting library 
	 		if(_currentinterface == 0)
	 		{
	 			tmpstr.push_back("Origin");
	 			_localstartinglibrary.push_back(tmpstr);		 			
	 		}
	 		else
	 		{
	 			tmpstr.push_back(_shootingconfigfile);
	 		}

	 		// Update index file of new configuration
	 		_indexcontents += tmpstr[0]+" "+tmpstr[1]+"\n";
	 		_indexfile << tmpstr[0] <<" "<<tmpstr[1]<<" "<<tmpstr[2]<<std::endl;
		}

		_currenthash++;
	}

	void ForwardFlux::ReadConfiguration(Snapshot* snapshot, std::string dumpfilename)
	{
		if(_world.rank() == 0)
			_dumpfilecontents = GetFileContents(dumpfilename.c_str());

		mpi::broadcast(_world, _dumpfilecontents, 0);

		auto& positions = snapshot->GetPositions();
		auto& velocities = snapshot->GetVelocities();
		auto& atomID = snapshot->GetAtomIDs();
		auto& forces = snapshot->GetForces();

		//Extract dump file information
		std::istringstream f(_dumpfilecontents);
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

			_shootingconfigfile = dumpfilename;
		}
	}

	void ForwardFlux::ClearFiles()
	{
		for(auto& tstr : _startinglibrary)
		{
			for(auto& tstr2 : tstr)
				tstr2.clear();
			tstr.clear();
		}

		_successes.clear();

		// close files and reopen to overwrite
		_indexfile.close();
		_indexfile.open(_indexfilename.c_str());
		_indexcontents = "";
			
		_resultsfile.close();
		_resultsfile.open(_resultsfilename.c_str());
		_resultscontents = "";

		_dumpfilecontents = "";

		_currentinterface = 0;
		_currenthash = 0;

		_fluxout = _fluxin = 0;

		_newrun = true;
	}

	// Check if starting a new run, if so generate new library 
	void ForwardFlux::SetUpNewRun(Snapshot* snapshot, const CVList& cvs)
	{
		int interface = AtInterface(cvs);
		// Get CV values check if at next interface, if so store configuration
		if(interface == 1 && _currentinterface == 0)
		{
			_currentinterface = interface;
			WriteConfiguration(snapshot);
			_successes[_currentinterface]++;
			_fluxout++;
		}
		else if(interface == 0 && _currentinterface == 1)
		{
			_currentinterface = interface;
			_fluxin++;
		}
		
		// See if found the required starting configs
		if(_localstartinglibrary.size() >= _requiredconfigs)
		{
			MPI_Barrier(_world);
			mpi::all_gather(_world, &_localstartinglibrary, _localstartinglibrary.size(), &_startinglibrary);
			
			_newrun = false;
			_restartfromlibrary = true;

			//Start at first location
			_currentstartingpoint = 0;
		}
	}

	// Pick a random configuration
	std::string ForwardFlux::PickConfiguration(int interface)
	{
		std::string configfilename;
		if(_world.rank() == 0)
		{
			std::vector<std::vector<std::string> > files; 
			if(!(ExtractInterfaceIndices(interface, files)))
			{
				std::cout<< "Could not locate any files at interface ";
				std::cout<< interface<<"!"<<std::endl;
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
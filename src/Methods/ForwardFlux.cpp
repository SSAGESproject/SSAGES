#include "ForwardFlux.h"

#include <iostream>

namespace mpi = boost::mpi;
namespace SSAGES
{
	void ForwardFlux::PreSimulation(Snapshot*, const CVList& cvs)
	{
		char file[1024];
				
		// Before the simulation begins:
		// Open up global file for indexing successful configurations
		if(_world.rank() == 0)
		{
			if(_NewRun)
			{
				sprintf(file, _indexfilename);
		 		_indexfile.open(file);

	 			sprintf(file, _resultsfilename);
	 		 	_resultsfile.open(file);
			}
			else
			{
				_indexcontents = GetFileContents(_indexfilename.c_str());
				_resultscontents = GetFileContents(_resultsfilename.c_str());

				sprintf(file, _indexfilename);
				_indexfile.open(file, std::ofstream::out | std::ofstream::app)

				sprintf(file, _resultsfilename);
				_resultsfile.open(file, std::ofstream::out | std::ofstream::app)

		 	}
		}

		if(!_NewRun)
		{
			mpi::broadcast(_world, _indexcontents, 0);
			mpi::broadcast(_world, _resultscontents, 0);
			mpi::broadcast(_comm, _dumpfilecontents, 0);
		}

		if(!_NewRun)
		{
			if(ExtractInterfaceIndices(_currentinterface, _startinglibrary))
			{
				if(_world.rank() == 0)
				{
					std::cout << "Could not locate an interfaces at " << _currentinterface << " in ";
					std::cout << _indexfilename << "! New starting configurations will be generated." << std::endl;
				}

				_currentinterface = 0;
				_NewRun = true;

				_indexcontents = "";
				for(auto& tstr : _startinglibrary)
					for(auto& tstr2 : tstr)
						tstr2.clear();
					tstr.clear();
			}
		}

		// Resize successes
		_successes.resize(_centers.size());
	}

	void ForwardFlux::PostIntegration(Snapshot* snapshot, const CVList& cvs)
	{
		//Check if starting a new run, if so generate new library 
		if(_NewRun)
		{
			int interface = AtInterface(cvs);
			// Get CV values check if at next interface, if so store configuration
			if((interface == 1 && _currentinterface == 0) || (interface == 0 && _currentinterface == 1))
			{
				_currenthash++;
				WriteConfiguration(snapshot);
				success[_currentinterface]++;
				_currentinterface = interface;
				_flux += interface;
			}

			// See if found the required starting configs
			if(_localstartinglibrary.size() >= _requiredconfigs)
			{
				MPI_Barrier(_world);
				mpi::all_gather(_world, _localstartinglibrary, _startinglibrary);
				_NewRun = false;
			}
		}

		
		// Check if at interface
		// If so:
		//		If success:
		// 			Write to dump file, store local interface number and hash
		//			increment success
		//		else:
		//			incrememnt failure
		//		wait for everyone else to finish
		// If not:
		//		Continue running
		//
		// If everyone is complete:
		//		From every walker update global indexing file
		//		Increment the current interface
		//		

		
		// If found a path to B stop simulation
		// NEED BETTER WAY OF DOING THIS
		if(ended)
		{
			PostSimulation(snapshot, cvs);
			_world.abort(0);
		}


	}

	void ForwardFlux::PostSimulation(Snapshot*, const CVList&)
	{
		// Calculate probabilites for each interface.
		// Write probabilities and paths to file

		//Close local and global files
		if(_world.rank() ==0)
		{
			_indexfile.close();
			_resultsfile.close();
		}
	}

}


/*
File Formats:
_indexfile
interface(some integer) dump_file_name(a string) hashID(a string)
example: 1 dump1.xyz 1-113




*/
#pragma once 

#include "Method.h"
#include "../CVs/CollectiveVariable.h"
#include <fstream>

namespace SSAGES
{
	// ForwardFlux sampling method 
	class ForwardFlux : public Method
	{
	private:

		// Output index file for storing information on where everything is.
		std::string _indexfilename;
		std::ofstream _indexfile;
		std::string _indexcontents;
		std::vector<std::vector<std::string> > _startinglibrary;
		std::vector<std::vector<std::string> > _localstartinglibrary;
		std::vector<std::vector<std::string> > _indexinformation;

		// Results file for end of simulation.
		std::string _resultsfilename;
		std::ofstream _resultsfile;
		std::string _resultscontents;

		// Location of the nodes to be used in determining FF interfaces
		std::vector<std::vector<double> > _centers;

		// User defined if we need to create a library of new configs or not
		bool _NewRun;

		// Current interface FF is shooting from
		int _currentinterface;

		// User defined number of starting configs needed per walker before starting FF
		int _requiredconfigs;

		// Number that keeps track of configs this walker has generated
		int _currenthash;

		// Name of file of configuration where shooting from
		std::string _shootingconfigfile;

		// Flux
		int _flux;

	public:
		// Create instance of Forward Flux with centers "centers". 
		ForwardFlux(boost::mpi::communicator& world,
				 boost::mpi::communicator& comm,
				 const std::vector<double>& centers,
				 unsigned int frequency) : 
		Method(frequency, world, comm), _centers(centers)
		{
		}

		// Pre-simulation hook.
		void PreSimulation(Snapshot* snapshot, const CVList& cvs) override;

		// Post-integration hook.
		void PostIntegration(Snapshot* snapshot, const CVList& cvs) override;

		// Post-simulation hook.
		void PostSimulation(Snapshot* snapshot, const CVList& cvs) override;
		
		void Serialize(Json::Value& json) const override
		{

		}

		// Extract all indices for a given interface. 
		// Return true if couldnt locate anything at a given interface
		bool ExtractInterfaceIndices(int interface, std::vector<std::vector<std::string> >& InterfaceIndices);
		{
			//Extract configuration indices for a given interface int
			std::istringstream f(_indexcontents);
			std::string line;
			while (std::getline(f, line))
			{
				string buf; // Have a buffer string
				stringstream ss(line); // Insert the string into a stream
				vector<string> tokens; // Create vector to hold our words

				while (ss >> buf)
				    tokens.push_back(buf);

				if(std::stoi(tokens[0]) == interface)
					InterfaceIndices.push_back(tokens);
			}

			if(InterfaceIndices.size() == 0)
				return true;

			return false;
		}

		int AtInterface(const CVList& cvs)
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

		void WriteConfiguration(Snapshot* snapshot)
		{
			const auto& positions = snapshot.GetPositions();
			const auto& velocities = snapshot.GetVelocities();
			const auto& atomID = snapshot.GetAtomIDs();

			if(_comm.rank() == 0)
			{	
	
				// Update index file of new configuration
				_indexcontents

				_shootingconfigfile
				// Write the dump file out
				std::string dumpfilename = "dump_"+std::to_string(_currentinterface)+"_"+std::to_string(_currenthash)+".dump";
				std::ofstream dumpfile;
				sprintf(file, dumpfilename);
		 		dumpfile.open(file);

		 		for(size_t i = 0; i< atomID.size(); i++)
		 		{
		 			dumpfile<<atiomID[i]<<" ";
		 			dumpfile<<positions[i][0]<<" "<<positions[i][1]<<" "<<positions[i][2]<<" "<<std::endl;
		 			dumpfile<<velocities[i][0]<<" "<<velocities[i][1]<<" "<<velocities[i][2]<<std::endl;
				}

		 		// Update starting library 
		 		if(_currentinterface = 0)
		 		{
		 			std::vector<std::string> tmpstr;
		 			tmpstr.push_back(std::to_string(_currentinterface));
		 			tmpstr.push_back(dumpfilename);
		 			tmpstr.push_back("Origin");

		 			_localstartinglibrary.push_back(tmpstr);
		 		}
			}
		}
	};
}
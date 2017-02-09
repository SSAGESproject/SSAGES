#pragma once 

#include "Method.h"
#include "../CVs/CollectiveVariable.h"
#include <fstream>

namespace SSAGES
{
	// Umbrella sampling method to constrain an arbitrary 
	// number of CVs at specified equilibrium distances.
	class Umbrella : public Method
	{
	private:
		// Vector of spring constants.
		std::vector<double> kspring_;

		// Vector of equilibrium distances.
		std::vector<double> centers_;

		// iterator for this method
		int currentiter_;

		// Output stream for umbrella data.
		std::ofstream umbrella_;

	public:
		// Create instance of umbrella with spring constants "kspring", 
		// and centers "centers". Note the sizes of the vectors should be 
		// commensurate with the number of CVs.
		Umbrella(boost::mpi::communicator& world,
				 boost::mpi::communicator& comm,
				 const std::vector<double>& kspring,
				 const std::vector<double>& centers,
				 unsigned int frequency) : 
		Method(frequency, world, comm), kspring_(kspring), centers_(centers)
		{}

		// Pre-simulation hook.
		void PreSimulation(Snapshot* snapshot, const CVList& cvs) override;

		// Post-integration hook.
		void PostIntegration(Snapshot* snapshot, const CVList& cvs) override;

		// Post-simulation hook.
		void PostSimulation(Snapshot* snapshot, const CVList& cvs) override;
		
		void PrintUmbrella(const CVList& cvs);

		void Serialize(Json::Value& json) const override
		{

		}

	};
}
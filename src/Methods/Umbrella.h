#pragma once 

#include "Method.h"
#include "../CVs/CollectiveVariable.h"

namespace SSAGES
{
	// Umbrella sampling method to constrain an arbitrary 
	// number of CVs at specified equilibrium distances.
	class Umbrella : public Method
	{
	private:
		// Vector of spring constants.
		std::vector<double> _kspring;

		// Vector of equilibrium distances.
		std::vector<double> _centers;

	public:
		// Create instance of umbrella with spring constants "kspring", 
		// and centers "centers". Note the sizes of the vectors should be 
		// commensurate with the number of CVs.
		Umbrella(boost::mpi::communicator& world,
				 boost::mpi::communicator& comm,
				 const std::vector<double>& kspring,
				 const std::vector<double>& centers,
				 unsigned int frequency) : 
		Method(frequency, world, comm), _kspring(kspring), _centers(centers)
		{}

		// Pre-simulation hook.
		void PreSimulation(Snapshot* snapshot, const CVList& cvs) override;

		// Post-integration hook.
		void PostIntegration(Snapshot* snapshot, const CVList& cvs) override;

		// Post-simulation hook.
		void PostSimulation(Snapshot* snapshot, const CVList& cvs) override;
	};
}
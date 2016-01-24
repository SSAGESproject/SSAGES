#pragma once 

#include "Method.h"
#include <iostream>
#include <iomanip>
#include "mpi.h"

namespace SSAGES
{
	class MockMethod : public Method
	{
	public:
		MockMethod(unsigned int frequency) : 
		Method(frequency)
		{
		}

		void PreSimulation(Snapshot* snapshot, const CVList& cvs) override
		{
		}

		void PostIntegration(Snapshot* snapshot, const CVList& cvs) override
		{
			using std::setw;
			using std::right;
			using std::setprecision;
			
			// An example of acquiring and printing some data from the snapshot.
			std::cout 
			<< setw(8) << right << snapshot->GetIteration() << std::fixed
			<< setw(13) << right << setprecision(4) << snapshot->GetVolume() 
			<< setw(13) << right << setprecision(7) << snapshot->GetPressure() 
			<< setw(13) << right << setprecision(7) << snapshot->GetEnergy()
			<< " // SSAGES mock method"
			<< std::endl; 
		}

		void PostSimulation(Snapshot* snapshot, const CVList& cvs) override
		{
			/*for(auto& a : snapshot->GetAtomIDs())
				std::cout << "Atom ID: " << a << std::endl;*/
		}
	};
}
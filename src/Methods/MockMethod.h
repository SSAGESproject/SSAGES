#pragma once 

#include "Method.h"
#include <iostream>
#include <iomanip>
#include <boost/mpi.hpp>
#include <fstream>

namespace SSAGES
{
	class MockMethod : public Method
	{
	private:
		std::ofstream _myout;

	public:
		MockMethod(boost::mpi::communicator& world,
				   boost::mpi::communicator& comm,
				   unsigned int frequency) : 
		Method(frequency,world,comm)
		{
			_myout.open("foo.out");
		}

		void PreSimulation(Snapshot*, const CVList&) override
		{
		}

		void PostIntegration(Snapshot* snapshot, const CVList& cvs) override
		{
			using std::setw;
			using std::right;
			using std::setprecision;
			

			std::cout << snapshot->GetWalkerID() 
			<< " " << snapshot->GetCommunicator().rank() << std::endl;
			auto& v = snapshot->GetForces();
			v[0][0] = 1.000000*v[0][0];
			/*
			// An example of acquiring and printing some data from the snapshot.
			std::cout 
			<< setw(8) << right << snapshot->GetIteration() << std::fixed
			<< setw(13) << right << setprecision(4) << snapshot->GetVolume() 
			<< setw(13) << right << setprecision(7) << snapshot->GetPressure() 
			<< setw(13) << right << setprecision(7) << snapshot->GetEnergy()
			<< " // SSAGES mock method"
			<< std::endl; */

			for(auto& cv : cvs)
				std::cout << cv->GetValue() << std::endl;
		}

		void PostSimulation(Snapshot*, const CVList&) override
		{
		}

		~MockMethod() { _myout.close(); }
	};
}
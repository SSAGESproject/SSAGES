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
		MockMethod(unsigned int frequency) : 
		Method(frequency)
		{
			_myout.open("foo.out");
		}

		void PreSimulation(Snapshot* snapshot, const CVList& cvs) override
		{
		}

		void PostIntegration(Snapshot* snapshot, const CVList& cvs) override
		{
			using std::setw;
			using std::right;
			using std::setprecision;
			
			boost::mpi::communicator comm;
			if(comm.rank() == 0)
			{
				for(auto& p : snapshot->GetPositions())
					_myout << p[0] << " " << p[1] << " " << p[2] << std::endl;

				_myout << std::endl;
 			}

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

		void PostSimulation(Snapshot* snapshot, const CVList& cvs) override
		{
		}

		~MockMethod() { _myout.close(); }
	};
}
#pragma once 

#include "Method.h"
#include <iostream>
#include <iomanip>
#include <boost/mpi.hpp>
#include <fstream>

namespace SSAGES
{
	//! Mock Sampling
	/*!
	 * \ingroup Methods
	 */
	class MockMethod : public Method
	{
	private:
		std::ofstream _myout; //!< Output stream.

	public:
		//! Constructor
		/*!
		 * \param world MPI global communicator.
		 * \param comm MPI local communicator.
		 * \param frequency Frequency with which this method is invoked.
		 */
		MockMethod(boost::mpi::communicator& world,
				   boost::mpi::communicator& comm,
				   unsigned int frequency) : 
		Method(frequency,world,comm)
		{
			_myout.open("foo.out");
		}

		//! Pre-simulation hook.
		void PreSimulation(Snapshot*, const CVList&) override
		{
		}

		//! Post-integration hook.
		/*!
		 * \param snapshot Current simulation snapshot.
		 * \param cvs List of CVs.
		 */
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

		//! Post-simulation hook.
		void PostSimulation(Snapshot*, const CVList&) override
		{
		}

		//! Destructor.
		~MockMethod() { _myout.close(); }
	};
}
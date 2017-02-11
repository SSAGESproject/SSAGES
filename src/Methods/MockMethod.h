/**
 * This file is part of
 * SSAGES - Suite for Advanced Generalized Ensemble Simulations
 *
 * Copyright 2016 Hythem Sidky <hsidky@nd.edu>
 *                Ben Sikora <bsikora906@gmail.com>
 *
 * SSAGES is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * SSAGES is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with SSAGES.  If not, see <http://www.gnu.org/licenses/>.
 */
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
		std::ofstream myout_; //!< Output stream.

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
			myout_.open("foo.out");
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
		~MockMethod() { myout_.close(); }
	};
}
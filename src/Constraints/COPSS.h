/**
 * This file is part of
 * SSAGES - Suite for Advanced Generalized Ensemble Simulations
 *
 * Copyright 2016 Jiyuan Li <jyli@uchicago.edu>
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

#include "Constraint.h"
#include <iostream>
#include <iomanip>
#include <boost/mpi.hpp>
#include <fstream>

namespace SSAGES
{
	//! Wrapper class for COPSS simulations.
	class COPSS : public Constraint
	{
	private:
		//! Output stream for debug output.
		std::ofstream myout_;

	public:
		//! Constructor.
		/*!
		 * \param comm MPI global communicator.
		 * \param frequency Frequency with which this constraint is called.
		 */
		COPSS(boost::mpi::communicator& comm,
				   unsigned int frequency) : 
		Constraint(frequency, comm)
		{
			myout_.open("test.out");
		}

		//! Pre-Simulation Hook.
		void PreSimulation(Snapshot*, const CVList&) override
		{
			myout_<<"In PreSimulation"<<std::endl;
		}

		//! Post-Integration Hook.
		void PostIntegration(Snapshot* snapshot, const CVList&) override
		{

			using std::setw;
			using std::right;
			using std::setprecision;
			

			myout_ << snapshot->GetWalkerID() 
			<< " " << snapshot->GetCommunicator().rank() << std::endl;
			auto& f = snapshot->GetForces();
			f[0][0] = 1.000000*f[0][0];
		}

		//! Post-Simulation Hook.
		void PostSimulation(Snapshot*, const CVList&) override
		{
			myout_ <<" Post simulation"<<std::endl;
		}

		//! Serialize the class
		/*!
		 * \warning Serialization not yet implemented.
		 */
		void Serialize(Json::Value& /* json */) const override
		{

		}

		//! Destructor.
		~COPSS() { myout_.close(); }
	};
}

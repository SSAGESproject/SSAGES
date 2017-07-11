/**
 * This file is part of
 * SSAGES - Suite for Advanced Generalized Ensemble Simulations
 *
 * Copyright 2017 Hythem Sidky <hsidky@nd.edu>
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

#include <mpi.h>

// Forward declare.
namespace Json {
    class Value;
}

namespace SSAGES
{
	//! Simulation driver.
    class Driver
    {
	private:
		//! Resource handler. 
		class ResourceHandler* rh_;

		//! Qbox hook. 
		class  QboxHook* qbhook_;

		//! MD iterations.
		int iterations_;

		//! QM iterations.
		int qmiterations_;

	public:
		Driver(class ResourceHandler* rh, class QboxHook* qbhook, int iter, int qmiter) : 
		rh_(rh), qbhook_(qbhook), iterations_(iter), qmiterations_(qmiter)
		{}

        void Run();

        //! Build a new Driver from JSON. 
		/*!
		 * \param json JSON root node containing driver (and children) specifications. 
		 * \param world MPI communicator containing all processors. 
		 * \return Pointer to newly created driver. 
		 * 
		 * \note Object lifetime is caller's responsibility!
		 */
		static Driver* Build(const Json::Value& json, const MPI_Comm& world);

		~Driver();
    };
}
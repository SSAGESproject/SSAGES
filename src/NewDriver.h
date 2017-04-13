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

#include "JSON/Serializable.h"
#include <mpi.h>
#include <vector>

namespace SSAGES
{
	//! Primary abstract class that drives a simulation. 
	/*!
	 * 
	 * Driver is responsible for intialzing the major components 
	 * of a simulation including the engine, methods, and collective
	 * variables. It is also responsible for the lifetime of the objects 
	 * it creates. 
	 * 
	 * Each simulation engine must implement a concrete driver which 
	 * handles the implementation-specific parts of the initialization 
	 * routines.
	 */
	class NewDriver : public Serializable
	{
	protected: 
		//! MPI communicator containing all processors. 
		MPI_Comm world_; 

		//! MPI communicator containing processors for specific walker. 
		MPI_Comm comm_;

		//! Walker ID for specific driver. 
		uint walkerid_; 

		//! Snapshot of system state (pointer). 
		class Snapshot* snapshot_;

		//! Vector of advanced sampling methods. 
		std::vector<class Method*> methods_;

		//! Collective variable manager. 
		class CVManager* cvmanager_;

	public:
		//! Constructor. 
		/*!
		 * \param world MPI communicator containing all processors. 
		 * \param comm MPI communicator containing walker-specific processors. 
		 * \param walkerid ID of the walker for the current processor.
		 * \param methods Vector of pointers to methods. 
		 * \param CVManager Pointer to CV manager. 
		 * 
		 * \note Driver will be responsible for lifetime of methods and CV manager. 
		 */
		NewDriver(const MPI_Comm& world, const MPI_Comm& comm, uint walkerid,
		          const std::vector<class Method*>& methods, class CVManager* cvmanager);

		//! Build a new Driver from JSON. 
		/*!
		 * \param json JSON root node containing driver (and children) specifications. 
		 * \param world MPI communicator containing all processors. 
		 * \return Pointer to newly created driver. 
		 * 
		 * \note Object lifetime is caller's responsibility!
		 */
		static NewDriver* Build(const Json::Value& json, const MPI_Comm& world);

        //! \copydoc Serializable::Serialize()
		void Serialize(Json::Value& json) const override;

		//! Destructor.
		~NewDriver();
	};
}
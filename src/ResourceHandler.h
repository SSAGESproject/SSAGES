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

#include <mxx/comm.hpp>
#include <vector>

// Forward declare.
namespace Json {
	class Value;
}

namespace SSAGES
{

	//! Class that handles SSAGES resources for a simulation.
	/*!
	 * 
	 * ResourceHandler is responsible for intialzing the major components 
	 * of a simulation including the engine, methods, and collective
	 * variables. It is also responsible for the lifetime of the objects 
	 * it creates. 
	 * 
	 * Each simulation engine must implement a driver which calls this 
	 * resource handler and passes it the appropriate hook. 
	 */
	class ResourceHandler
	{
	private: 
		//! MPI communicator containing all processors. 
		mxx::comm world_; 

		//! MPI communicator containing processors for specific walker. 
		mxx::comm comm_;

		//! Walker ID for specific driver. 
		uint walkerid_; 

		//! Number of walkers.
		uint nwalkers_;

		//! Snapshot of system state (pointer). 
		class Snapshot* snapshot_;

		//! Vector of advanced sampling methods. 
		std::vector<class Method*> methods_;

		//! Collective variable manager. 
		class CVManager* cvmanager_;

		//! Hook pointer.
		class Hook* hook_;

		//! Input file vector. 
		std::vector<std::string> inputs_;

	public:
		//! Constructor. 
		/*!
		 * \param world MPI communicator containing all processors. 
		 * \param comm MPI communicator containing walker-specific processors. 
		 * \param walkerid ID of the walker for the current processor.
		 * \param methods Vector of pointers to methods. 
		 * \param CVManager Pointer to CV manager. 
		 * 
		 * \note ResourceHandler will be responsible for lifetime of methods and CV manager. 
		 */
		ResourceHandler(mxx::comm&& world, mxx::comm&& comm, uint walkerid,
		          const std::vector<class Method*>& methods, class CVManager* cvmanager);

		
		std::string GetInput() const
		{
			return inputs_[walkerid_];
		}

		MPI_Comm GetLocalComm() const
		{
			return comm_;
		}

		MPI_Comm GetWorldComm() const
		{
			return world_;
		}

		const mxx::comm& GetLocalMxxComm() const
		{
			return comm_;
		}
		
		const mxx::comm& GetWorldMxxComm() const
		{
			return world_;
		}

		uint GetWalkerID() const
		{
			return walkerid_;
		}

		uint GetNumWalkers() const 
		{
			return nwalkers_;
		}

		void ConfigureHook(class Hook* hook);

		//! Build a new ResourceHandler from JSON. 
		/*!
		 * \param json JSON root node containing driver (and children) specifications. 
		 * \param world MPI communicator containing all processors. 
		 * \return Pointer to newly created driver. 
		 * 
		 * \note Object lifetime is caller's responsibility!
		 */
		static ResourceHandler* Build(const Json::Value& json, const MPI_Comm& world);

		//! Destructor.
		~ResourceHandler();
	};
}
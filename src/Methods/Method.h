/**
 * This file is part of
 * SSAGES - Suite for Advanced Generalized Ensemble Simulations
 *
 * Copyright 2016 Ben Sikora <bsikora906@gmail.com>
 *                Hythem Sidky <hsidky@nd.edu>
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

#include "../EventListener.h"
#include "../JSON/Serializable.h"
#include <mxx/comm.hpp>
#include <functional>

// Forward declare.
namespace Json {
	class Value;
}

namespace SSAGES
{
	//! Interface for Method implementations.
	/*!
	 * The base method class from which advanced sampling routines derive. 
	 * A method is allowed to manipulate a simulation at three points: 
	 * before the simulation begins (usually initialization), after each 
	 * integration step by the simulation engine, and after the integration 
	 * steps are complete (usually cleanup). 
	 *
	 * \ingroup Methods
	 */
	class Method : public EventListener, public Serializable
	{
	protected:
		mxx::comm world_; //!< Global MPI communicator
		mxx::comm comm_; //!< Local MPI communicator

	public:
		//! Constructor
		/*!
		 * \param frequency Frequency of sampling.
		 * \param world Global MPI communicator.
		 * \param comm MPI communicator of walker.
		 *
		 * Frequency of sampling must be specified by all methods.
		 */
		Method(uint frequency, const MPI_Comm& world, const MPI_Comm& comm) : 
		EventListener(frequency), world_(world), comm_(comm)
		{}

		//! Method call prior to simulation initiation.
		/*!
		 * \param snapshot Pointer to the simulation snapshot.
		 * \param cvs List of CVs.
		 *
		 * This function will be called before the simulation is started.
		 */
		virtual void PreSimulation(Snapshot* snapshot, const CVList& cvs) override = 0;

		//! Method call post integration.
		/*!
		 * \param snapshot Pointer to the simulation snapshot.
		 * \param cvs List of CVs.
		 *
		 * This function will be called after each integration step.
		 */
		virtual void PostIntegration(Snapshot* snapshot, const CVList& cvs) override = 0;

		//! Method call post simulation.
		/*!
		 * \param snapshot Pointer to the simulation snapshot.
		 * \param cvs List of CVs.
		 *
		 * This function will be called after the end of the simulation run.
		 */
		virtual void PostSimulation(Snapshot* snapshot, const CVList& cvs) override = 0;
		
		//! Build a derived method from JSON node.
		/*!
		 * \param json JSON Value containing all input information.
		 * \param world MPI global communicator.
		 * \param comm MPI local communicator.
		 * \param path Path for JSON path specification.
		 * \return Pointer to the Method built. nullptr if an unknown error occurred.
		 *
		 * This function builds a registered method from a JSON node. The difference
		 * between this function and "Build" is that this automatically determines the 
		 * appropriate derived type based on the JSON node information.
		 *
		 * \note Object lifetime is the caller's responsibility.
		 */
		static Method* BuildMethod(const Json::Value& json, 
		                           const MPI_Comm& world, 
		                           const MPI_Comm& comm, 
		                           const std::string& path);

		//! Destructor
		virtual ~Method() 
		{
		}
	};
}

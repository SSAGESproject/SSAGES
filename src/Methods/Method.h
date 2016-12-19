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
#include "../Grid.h"
#include <boost/mpi.hpp>

// Forward declare.
namespace Json {
	class Value;
}

namespace SSAGES
{
	//! Interface for Method implementations.
	/*!
	 * This class is designed to be the base class on which all Metadynamics
	 * methods are based.
	 *
	 * \ingroup Methods
	 */
	class Method : public EventListener, public Serializable
	{
	protected:
		boost::mpi::communicator world_; //!< Global MPI communicator
		boost::mpi::communicator comm_; //!< Local MPI communicator

		//! Pointer to grid
		Grid<int>* grid_;

		//! Number of the method iteration.
		unsigned int iteration_;

	public:
		//! Constructor
		/*!
		 * \param frequency Frequency of sampling.
		 * \param world MPI world communicator.
		 * \param comm MPI local communicator.
		 *
		 * Frequency of sampling must be specified by all methods.
		 */
		Method(unsigned int frequency, 
			boost::mpi::communicator& world, 
			boost::mpi::communicator& comm) : 
		EventListener(frequency), world_(world), comm_(comm),
		grid_(nullptr), iteration_(0){}

		//! Set Method's iteration.
		/*!
		 * \param iter int value for what the method iteration should be.
		 */
		void SetIteration(int iter) {iteration_ = iter;}

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

		//! Set up the Grid
		/*!
		 * \param json JSON Value containing the input information
		 * \param path Path for JSON path specification.
		 */
		void BuildGrid(const Json::Value& json, const std::string& path)
		{
			grid_ = Grid<int>::BuildGrid(json, path);
		}

		//! Get the Grid
		Grid<int>* GetGrid() const {return grid_;}

		//! Set up the Method
		/*!
		 * \param json JSON Value containing all input information.
		 * \param world MPI global communicator.
		 * \param comm MPI local communicator.
		 * \param path Path for JSON path specification.
		 * \return Pointer to the Method built. nullptr if an unknown error occurred.
		 *
		 * This function builds a method from a JSON node. It will return a
		 * nullptr when an unknown error occurred, but generally it will throw
		 * a BuildException on failure.
		 *
		 * \note Object lifetime is the caller's responsibility.
		 */
		static Method* BuildMethod(const Json::Value& json,
								boost::mpi::communicator& world, 
								boost::mpi::communicator& comm,
							   	const std::string& path);

		//! Destructor
		virtual ~Method() 
		{
			delete grid_;
		}
	};
}

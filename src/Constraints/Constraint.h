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

#include "../EventListener.h"
#include "../Snapshot.h"
#include <boost/mpi.hpp>

// Forward declare.
namespace Json {
	class Value;
}

namespace SSAGES
{
	// Forward declare.
	class Constraint;

	//! List of Constraints
	using ConstraintList = std::vector<Constraint*>;

	//! Interface for Constraint implementations.
	class Constraint : public EventListener, public Serializable
	{
	protected:
		//! MPI global communicator.
		boost::mpi::communicator comm_;

	public:
		//! Constructor
		/*!
		 * \param frequency Frequency of sampling.
		 * \param comm MPI global communicator.
		 *
		 * \note Frequency of sampling must be specified by all methods.
		 */
		Constraint(unsigned int frequency,  
			boost::mpi::communicator& comm) : 
		EventListener(frequency), comm_(comm){}

		//! Destructor
		virtual ~Constraint(){}

		//! Method call prior to simulation initiation.
		/*!
		 * \param snapshot Simulation snapshot.
		 * \param cvs List of CVs.
		 */
		virtual void PreSimulation(Snapshot* snapshot, const CVList& cvs) override = 0;

		//! Method call post integration.
		/*!
		 * \param snapshot Pointer to current snapshot.
		 * \param cvs List of CVs.
		 */
		virtual void PostIntegration(Snapshot* snapshot, const CVList& cvs) override = 0;

		//! Method call post simulation.
		/*!
		 * \param snapshot Pointer to current snapshot.
		 * \param cvs List of CVs.
		 */
		virtual void PostSimulation(Snapshot* snapshot, const CVList& cvs) override = 0;

		//! Build a constraint from a JSON node.
		/*!
		 * \param json JSON value containing input information.
		 * \param comm MPI communicator.
		 * \return Pointer to the new constraint. \c nullptr in case of an
		 *         unknown error.
		 *
		 * This function builds a constraint from a JSON node. Returns a pointer
		 * to the built constraint. If an unknown error occured, the return value
		 * is \c nullptr, but in general, a BuildException will be thrown on
		 * failure.
		 *
		 * \note Object lifetime is the caller's responsibility.
		 */
		static Constraint* BuildConstraint(const Json::Value& json, 
							boost::mpi::communicator& comm);

		//! Overloaded function allowing JSON path specification.
		/*!
		 * \param json JSON value containing input information.
		 * \param comm MPI global communicator.
		 * \param path Path for JSON path specification.
		 * \return Pointer to the new constraint. \c nullptr in case of an
		 *         unknown error.
		 */
		static Constraint* BuildConstraint(const Json::Value& json,
							boost::mpi::communicator& comm, 
							   const std::string& path);

		//! Build constraint.
		/*!
		 * \param json JSON value containing input information.
		 * \param clist List of constraints.
		 * \param comm MPI global communicator.
		 * \param path Path for JSON path specification.
		 *
		 * This function builds a new constraint and adds it to the specified
		 * list of constraints. On failure, an exception is thrown.
		 * \note Object lifetime management is caller's responsibility.
		 */
		static void BuildConstraint(const Json::Value& json, 
							   ConstraintList& clist,
							   boost::mpi::communicator& comm,
							   const std::string& path);
	};
}

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

	// Typedefs
	using ConstraintList = std::vector<Constraint*>;

	// Interface for Constraint implementations.
	class Constraint : public EventListener, public Serializable
	{
	protected:
		boost::mpi::communicator _comm;

	public:
		// Frequency of sampling must be specified by all methods.
		Constraint(unsigned int frequency,  
			boost::mpi::communicator& comm) : 
		EventListener(frequency), _comm(comm){}

		virtual ~Constraint(){}

		// Method call prior to simulation initiation.
		virtual void PreSimulation(Snapshot* snapshot, const CVList& cvs) override = 0;

		// Method call post integration.
		virtual void PostIntegration(Snapshot* snapshot, const CVList& cvs) override = 0;

		// Method call post simulation.
		virtual void PostSimulation(Snapshot* snapshot, const CVList& cvs) override = 0;

		// Builds a constraint from a JSON node. Returns a pointer to the built constraint.
		// If return value is nullptr, 
		// then an unknown error occurred. It will throw a BuildException on failure. 
		// Object lifetime is the caller's responsibility. 
		static Constraint* BuildConstraint(const Json::Value& json, 
							boost::mpi::communicator& comm);

		// Overloaded function allowing JSON path specification.
		static Constraint* BuildConstraint(const Json::Value& json,
							boost::mpi::communicator& comm, 
							   const std::string& path);

		// Builds constraints and adds them to the constraint list. 
		// Throws exception on failure. 
		// Object lifetime management is caller's responsibility. 
		static void BuildConstraint(const Json::Value& json, 
							   ConstraintList& clist,
							   boost::mpi::communicator& comm,
							   const std::string& path);
	};
}
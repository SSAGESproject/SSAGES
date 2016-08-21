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
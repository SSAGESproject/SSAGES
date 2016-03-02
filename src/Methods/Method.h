#pragma once

#include "../EventListener.h"
#include <boost/mpi.hpp>

namespace SSAGES
{
	// Interface for Method implementations.
	class Method : public EventListener
	{
	protected:
		boost::mpi::communicator _world, _comm;

	public:
		// Frequency of sampling must be specified by all methods.
		Method(unsigned int frequency, 
			boost::mpi::communicator& world, 
			boost::mpi::communicator& comm) : 
		EventListener(frequency), _world(world), _comm(comm) {}

		// Method call prior to simulation initiation.
		virtual void PreSimulation(Snapshot* snapshot, const CVList& cvs) override = 0;

		// Method call post integration.
		virtual void PostIntegration(Snapshot* snapshot, const CVList& cvs) override = 0;

		// Method call post simulation.
		virtual void PostSimulation(Snapshot* snapshot, const CVList& cvs) override = 0;

		virtual ~Method() {}
	};
}
#pragma once

#include "../JSON/Serializable.h"
#include "json/json.h"
#include <vector>
#include <boost/mpi.hpp>

namespace SSAGES
{
	// Base class for standard simulation. An simulation is provided with a 
	// reference to a World and a ForceFieldManager. The World is responsible 
	// for handling the "box" geometry and particles. The ForeceFieldManager 
	// contains the forcefield data for all Particle types and interactions.
	class Driver
	{
	private:

		CVList _cvs;
		Simulation* _sim;
		Method* _method;
		mpi::communicator _world_comm;
		mpi::communicator _local_comm;

	public:

		virtual void Run() = 0;

		virtual void Build();

		// Serialize
		virtual void Serialize(Json::Value& json) const override
		{

		}

		/* Static Builder methods */

		// Builds a simulation from a JSON node, along with the
		// required dependencies. If a dependency is not needed/initialized, 
		// simply pass a nullptr. I
		static Simulation* BuildSimulation(const Json::Value& json,
										   WorldManager* wm, 
										   ForceFieldManager* ffm,
										   MoveManager* mm, 
										   DOSOrderParameter* dop,
										   Histogram* hist);
		virtual ~Simulation(){}
	};
}

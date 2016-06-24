#pragma once

#include "../EventListener.h"
#include "../Grids/Grid.h"
#include <boost/mpi.hpp>

// Forward declare.
namespace Json {
	class Value;
}

namespace SSAGES
{
	// Interface for Method implementations.
	class Method : public EventListener, public Serializable
	{
	protected:
		boost::mpi::communicator _world, _comm;

		Grid* _grid;

	public:
		// Frequency of sampling must be specified by all methods.
		Method(unsigned int frequency, 
			boost::mpi::communicator& world, 
			boost::mpi::communicator& comm) : 
		EventListener(frequency), _world(world), _comm(comm){}

		// Method call prior to simulation initiation.
		virtual void PreSimulation(Snapshot* snapshot, const CVList& cvs) override = 0;

		// Method call post integration.
		virtual void PostIntegration(Snapshot* snapshot, const CVList& cvs) override = 0;

		// Method call post simulation.
		virtual void PostSimulation(Snapshot* snapshot, const CVList& cvs) override = 0;

		void BuildGrid(const Json::Value& json, const std::string& path)
		{
			_grid = Grid::BuildGrid(json, path);
		}

		// Builds a method from a JSON node. Returns a pointer to the built Method.
		// If return value is nullptr, 
		// then an unknown error occurred. It will throw a BuildException on failure. 
		// Object lifetime is the caller's responsibility. 
		// Function allowing JSON path specification.
		static Method* BuildMethod(const Json::Value& json,
								boost::mpi::communicator& world, 
								boost::mpi::communicator& comm,
							   	const std::string& path);

		virtual ~Method() 
		{
			delete _grid;
		}
	};
}
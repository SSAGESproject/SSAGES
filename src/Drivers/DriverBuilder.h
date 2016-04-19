#pragma once 

#include "json/json.h"
#include "Simulation/SimException.h"
#include "CVs/CollectiveVariable.h"
#include "Methods/Method.h"
#include "SimObserver.h"
#include "Simulation.h"
#include <iostream>
#include <iomanip>

namespace SSAGES
{
	// Class for creating, managing and running a simulation. 
	class SimBuilder
	{
	private: 
		CVList _cvs;
		Simulation* _sim;
		Method* _method;
		mpi::communicator _world_comm;
		mpi::communicator _local_comm;




		int _ltot, _msgw, _notw;

	public:
		SimBuilder() : 
		_cvs(), _sim(nullptr), _method(nullptr), 
		_ltot(81), _msgw(51), _notw(_ltot - _msgw)
		{}
		
		bool BuildSimulation(const std::string& filename);
		Simulation* GetSimulation() { return _sim; }
		
		~SimBuilder()
		{
			for(auto& cv : _cvs)
				delete cv;
			_worlds.clear();
		}
	};
}
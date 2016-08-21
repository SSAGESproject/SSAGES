#pragma once 

#include "Constraint.h"
#include <iostream>
#include <iomanip>
#include <boost/mpi.hpp>
#include <fstream>

namespace SSAGES
{
	class COPSS : public Constraint
	{
	private:
		std::ofstream _myout;

	public:
		COPSS(boost::mpi::communicator& comm,
				   unsigned int frequency) : 
		Constraint(frequency, comm)
		{
			_myout.open("test.out");
		}

		void PreSimulation(Snapshot*, const CVList&) override
		{
			_myout<<"In PreSimulation"<<std::endl;
		}

		void PostIntegration(Snapshot* snapshot, const CVList&) override
		{

			using std::setw;
			using std::right;
			using std::setprecision;
			

			_myout << snapshot->GetWalkerID() 
			<< " " << snapshot->GetCommunicator().rank() << std::endl;
			auto& f = snapshot->GetForces();
			f[0][0] = 1.000000*f[0][0];
		}

		void PostSimulation(Snapshot*, const CVList&) override
		{
			_myout <<" Post simulation"<<std::endl;
		}

		void Serialize(Json::Value& json) const override
		{

		}

		~COPSS() { _myout.close(); }
	};
}
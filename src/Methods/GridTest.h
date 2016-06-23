#pragma once

#include "Method.h"
#include "../CVs/CollectiveVariable.h"
#include <fstream>
#include <vector>


namespace SSAGES
{

	class GridTest : public Method
	{

	public: 
		GridTest(boost::mpi::communicator& world,
					boost::mpi::communicator& comm,
					unsigned int frequency): 
		Method(frequency, world, comm)
		{
		}

		// Pre-simulation hook.
		void PreSimulation(Snapshot* snapshot, const CVList& cvs) override
		{
			_grid->PrintGrid();
		}

		// Post-integration hook.
		void PostIntegration(Snapshot* snapshot, const CVList& cvs) override
		{

			std::vector<int> Indices;
			int NDim = _grid->GetDimension();
			for(int i=0; i<NDim; i++)
				Indices.push_back(2);

			float val = 10;
			val += _grid->GetValue(Indices);
			_grid->PrintGrid();
			_grid->SetValue(Indices, val);
			_grid->PrintGrid();

			std::vector<double> spacing = _grid->GetSpacing();
			std::vector<double> upper = _grid->GetUpper();
			std::vector<double> lower = _grid->GetLower();
			std::vector<bool> periodic = _grid->GetPeriodic();
			std::vector<int> numpoints = _grid->GetNumPoints();

			for(int i = 0;i<NDim;i++)
				std::cout<<spacing[i]<<" "<<upper[i]<<" "<<lower[i]<<" "<<periodic[i]<<" "<<numpoints[i]<<std::endl;

			std::cout<<spacing.size()<<" "<<upper.size()<<" "<<lower.size()<<" "<<periodic.size()<<" "<<numpoints.size()<<std::endl;
			
		}

		// Post-simulation hook.
		void PostSimulation(Snapshot* snapshot, const CVList& cvs) override
		{

		}

		void Serialize(Json::Value& json) const override
		{

		}

		~GridTest() {}
	};
}
			

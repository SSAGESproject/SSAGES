/**
 * This file is part of
 * SSAGES - Suite for Advanced Generalized Ensemble Simulations
 *
 * Copyright 2016 Ben Sikora <bsikora906@gmail.com>
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

#include "Method.h"
#include "../CVs/CollectiveVariable.h"
#include <fstream>
#include <vector>


namespace SSAGES
{

	//! Test class for the Grid
	/*!
	 * \ingroup Methods
	 */
	class GridTest : public Method
	{

	public: 
		//! Constructor
		/*!
		 * \param world MPI global communicator.
		 * \param comm MPI local communicator.
		 * \param frequency Frequency with which this method is invoked.
		 */
		GridTest(boost::mpi::communicator& world,
					boost::mpi::communicator& comm,
					unsigned int frequency): 
		Method(frequency, world, comm)
		{
		}

		//! Pre-simulation hook.
		void PreSimulation(Snapshot* /* snapshot */, const CVList& /* cvs */) override
		{
			if(!grid_)
			{
				throw BuildException({"Method expected a grid but no grid built!"});
				world_.abort(-1);
			}
			grid_->PrintGrid();

			iteration_ = 0;
		}

		//! Post-integration hook.
		void PostIntegration(Snapshot* /* snapshot */, const CVList& /* cvs */) override
		{
			std::vector<int> Indices;
			std::vector<double> vals;
			int NDim = grid_->GetDimension();
			for(int i=0; i<NDim; i++)
			{
				Indices.push_back(2);
				vals.push_back(0);
			}

			double val = 10;
			val += grid_->GetValue(Indices);
			grid_->PrintGrid();
			grid_->SetValue(Indices, val);
			grid_->PrintGrid();

			std::vector<double> spacing = grid_->GetSpacing();
			std::vector<double> upper = grid_->GetUpper();
			std::vector<double> lower = grid_->GetLower();
			std::vector<bool> periodic = grid_->GetPeriodic();
			std::vector<int> numpoints = grid_->GetNumPoints();

			for(int i = 0;i<NDim;i++)
				std::cout<<spacing[i]<<" "<<upper[i]<<" "<<lower[i]<<" "<<periodic[i]<<" "<<numpoints[i]<<std::endl;

			std::cout<<spacing.size()<<" "<<upper.size()<<" "<<lower.size()<<" "<<periodic.size()<<" "<<numpoints.size()<<std::endl;
			std::cout<<"Indices at value 0 are: ";
			Indices = grid_->GetIndices(vals);
			for(int i = 0; i<NDim; i++)
				std::cout<<Indices[i]<<" ";
			std::cout<<std::endl;

			iteration_ ++;
		}

		//! Post-simulation hook.
		void PostSimulation(Snapshot* /* snapshot */, const CVList& /* cvs */) override
		{

		}

		//! \copydoc Serializable::Serialize()
		/*!
		 * \warning Serialization not implemented yet!
		 */
		void Serialize(Json::Value& json) const override
		{
			json["type"] = "GridTest";

		}

		//! Destructor.
		~GridTest() {}
	};
}
			

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
		/*!
		 * \param snapshot Current simulation snapshot.
		 * \param cvs List of CVs.
		 */
		void PreSimulation(Snapshot* snapshot, const CVList& cvs) override
		{
			if(!_grid)
			{
				throw BuildException({"Method expected a grid but no grid built!"});
				_world.abort(-1);
			}
			_grid->PrintGrid();

			_iteration = 0;
		}

		//! Post-integration hook.
		/*!
		 * \param snapshot Current simulation snapshot.
		 * \param cvs List of CVs.
		 */
		void PostIntegration(Snapshot* snapshot, const CVList& cvs) override
		{

			std::vector<int> Indices;
			std::vector<double> vals;
			int NDim = _grid->GetDimension();
			for(int i=0; i<NDim; i++)
			{
				Indices.push_back(2);
				vals.push_back(0);
			}

			double val = 10;
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
			std::cout<<"Indices at value 0 are: ";
			Indices = _grid->GetIndices(vals);
			for(int i = 0; i<NDim; i++)
				std::cout<<Indices[i]<<" ";
			std::cout<<std::endl;

			_iteration ++;
		}

		//! Post-simulation hook.
		/*!
		 * \param snapshot Current simulation snapshot.
		 * \param cvs List of CVs.
		 */
		void PostSimulation(Snapshot* snapshot, const CVList& cvs) override
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
			

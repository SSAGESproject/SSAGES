/**
 * This file is part of
 * SSAGES - Suite for Advanced Generalized Ensemble Simulations
 *
 Copyright 2016 Ben Sikora <bsikora906@gmail.com>
 *              Ashley Guo <azguo@uchicago.edu>
 *              Cody Bezik <bezik@uchicago.edu>
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

#include <numeric>

#include "Method.h"
#include <fstream>

namespace SSAGES
{
	//! String base class for FTS, Swarm, and elastic band
	/*!
	 * \ingroup Methods
	 *
	 * Implementation of a multi-walker finite string
	 * method with hard wall voronoi cells and running block averages.
	 */
	class StringMethod : public Method
	{
	protected:	
		
		//! CV starting location values.
		std::vector<double> centers_;

		//! CV starting location values.
		std::vector<double> newcenters_;

		//! The world's strings centers for each CV.
		/*!
		 * worldstring_[node#][cv#]
		 */
		std::vector<std::vector<double> > worldstring_;

		//! The node this belongs to
		int mpiid_;

		//! Tolerance criteria for determining when to stop string (default 0 if no tolerance criteria)
		std::vector<double> tol_;

		//! Number of nodes on a string
		int numnodes_;

		//! Maximum cap on number of string method iterations performed
		unsigned int maxiterator_;

		//! Vector of spring constants.
		std::vector<double> cvspring_;

		//! The local method iterator
		unsigned int iterator_;

		//! The global method iteration.
		uint iteration_;

		//! Output stream for string data.
		std::ofstream stringout_;

		//! Neighbor to send info to.
		int sendneigh_;

		//! Neighbor to gain info from.
		int recneigh_;

		//! Store positions for starting trajectories
		std::vector<std::vector<double>> prev_positions_;

		//! Store velocities for starting trajectories
		std::vector<std::vector<double>> prev_velocities_;

		//! Store atom IDs for starting trajectories
		std::vector<std::vector<int>> prev_IDs_;

		//! Updates the position of the string.
		virtual void StringUpdate() = 0;

		//! Helper function for calculating distances
		/*!
		 * \param x List of coordinates.
		 * \param y List of coordinates.
		 * \return Sum of distances between the x-values and y-values.
		 */
		double distance(const std::vector<double>& x, const std::vector<double>& y) const
		{
			double distance = 0;
			for (size_t i = 0; i < x.size(); i++)
				distance += pow((x[i] - y[i]),2);

			return sqrt(distance);
		}

		//! Prints the CV positions to file
		void PrintString(const CVList& CV);

		//! Gather neighbors over MPI
		/*!
		 * \param lcv0 Pointer to store array of lower CV values.
		 * \param ucv0 Pointer to store array of upper CV values.
		 */
		void GatherNeighbors(std::vector<double> *lcv0, std::vector<double> *ucv0);

		//! Reparameterize the string
		/*!
		 * \param alpha_star Factor for reparametrization.
		 */
		void StringReparam(double alpha_star);

		//! Update the world string over MPI
		/*!
		 * \param cvs List of CVs.
		 */
		void UpdateWorldString(const CVList& cvs);

		//! Check whether tolerance criteria has been met.
		bool TolCheck() const
		{
			if(tol_.size() == 0)
				return false;

			for(size_t i = 0; i < centers_.size(); i++)
            {
                if(fabs(centers_[i] - worldstring_[mpiid_][i]) > tol_[i])
                {
					return false;
                }
            }

			return true;
		}

		//! Check if method reached one of the exit criteria.
		/*!
		 * \param CV list of CVs.
		 * \return True if one of the exit criteria is met.
		 *
		 * The string method exits if either the maximum number of iteration has
		 * been reached or if all CVs are within the given tolerance thresholds.
		 */
		bool CheckEnd(const CVList& CV);

	public:
		//! Constructor
		/*!
		 * \param world MPI global communicator.
		 * \param comm MPI local communicator.
		 * \param centers List of centers.
		 * \param maxiterations Maximum number of iterations.
		 * \param cvspring Spring constants for cvs.
		 * \param frequency Frequency with which this method is invoked.
		 */
		StringMethod(const MPI_Comm& world,
					 const MPI_Comm& comm,
					 const std::vector<double>& centers,
					 unsigned int maxiterations,
					 const std::vector<double> cvspring,
			 		 unsigned int frequency) : 
		Method(frequency, world, comm), centers_(centers), 
		maxiterator_(maxiterations), 
		cvspring_(cvspring), iterator_(1), iteration_(0)
		{
			newcenters_.resize(centers_.size(), 0);
		}

		//! Pre-simulation hook.
		void PreSimulation(Snapshot* snapshot, const class CVManager& cvmanager) override;

		//! Post-integration hook.
		virtual void PostIntegration(Snapshot* snapshot, const class CVManager& cvmanager) override = 0;

		// Post-simulation hook.
		void PostSimulation(Snapshot*, const class CVManager&) override
		{
			stringout_.close();
		}

		//! Set the tolerance for quitting method
		/*!
		 * \param tol List of tolerances for each CV.
		 *
		 * Set the tolerance levels until which the method should run. The method
		 * quit, when the tolerance level is reached for all CVs.
		 */
		void SetTolerance(std::vector<double> tol)
		{
			for(auto& t : tol)
				tol_.push_back(t);
		}

		//! Communicate neighbor lists over MPI
        void SetSendRecvNeighbors();

		//! \copydoc Method::BuildMethod()
		static StringMethod* Build(const Json::Value& json, 
		                               const MPI_Comm& world,
		                               const MPI_Comm& comm,
					                   const std::string& path);

		//! Destructor
		virtual ~StringMethod() {}
	};
}
			

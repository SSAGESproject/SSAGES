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
#include "../CVs/CollectiveVariable.h"
#include <fstream>
#include "../spline.h"

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
	private:

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
		void PrintString(const CVList& CV)
		{
			if(comm_.rank() == 0)
			{
		        //Write node, iteration, centers of the string and current CV value to output file
		        stringout_.precision(8);
		        stringout_ << mpiid_ << " " << iteration_ << " ";

		        for(size_t i = 0; i < centers_.size(); i++)
		            stringout_ << worldstring_[mpiid_][i] << " " << CV[i]->GetValue() << " ";

		        stringout_ << std::endl;
		    }
		}

		//! Gather neighbors over MPI
		/*!
		 * \param lcv0 Pointer to store array of lower CV values.
		 * \param ucv0 Pointer to store array of upper CV values.
		 */
		void GatherNeighbors(std::vector<double> *lcv0, std::vector<double> *ucv0)
		{
			MPI_Status status;

			if(comm_.rank() == 0)
			{
				MPI_Sendrecv(&centers_[0], centers_.size(), MPI_DOUBLE, sendneigh_, 1234,
					&(*lcv0)[0], centers_.size(), MPI_DOUBLE, recneigh_, 1234, 
					world_, &status);

				MPI_Sendrecv(&centers_[0], centers_.size(), MPI_DOUBLE, recneigh_, 4321,
			       &(*ucv0)[0], centers_.size(), MPI_DOUBLE, sendneigh_, 4321, 
			       world_, &status);
			}

			MPI_Bcast(&(*lcv0)[0],centers_.size(),MPI_DOUBLE,0,comm_);
			MPI_Bcast(&(*ucv0)[0],centers_.size(),MPI_DOUBLE,0,comm_);
		}

		//! Reparameterize the string
		/*!
		 * \param alpha_star Factor for reparametrization.
		 */
		void StringReparam(double alpha_star)
		{
			std::vector<double> alpha_star_vector(numnodes_,0.0);

			//Reparameterization
			//Alpha star is the uneven mesh, approximated as linear distance between points
			if(comm_.rank()==0)
				alpha_star_vector[mpiid_] = mpiid_ == 0 ? 0 : alpha_star;

			//Gather each alpha_star into a vector 
			MPI_Allreduce(MPI_IN_PLACE, &alpha_star_vector[0], numnodes_, MPI_DOUBLE, MPI_SUM, world_);

			for(size_t i = 1; i < alpha_star_vector.size(); i++)
			    alpha_star_vector[i] += alpha_star_vector[i-1];
			
			for(size_t i = 1; i < alpha_star_vector.size(); i++)
			    alpha_star_vector[i] /= alpha_star_vector[numnodes_ - 1];

			tk::spline spl; //Cubic spline interpolation

			for(size_t i = 0; i < centers_.size(); i++)
			{
				std::vector<double> cvs_new(numnodes_, 0.0);

				if(comm_.rank() == 0)
					cvs_new[mpiid_] = centers_[i];

				MPI_Allreduce(MPI_IN_PLACE, &cvs_new[0], numnodes_, MPI_DOUBLE, MPI_SUM, world_);

			    spl.set_points(alpha_star_vector, cvs_new);
			    centers_[i] = spl(mpiid_/(numnodes_ - 1.0)); 
			}
		}

		//! Update the world string over MPI
		/*!
		 * \param cvs List of CVs.
		 */
		void UpdateWorldString(const CVList& cvs)
		{
			for(size_t i = 0; i < centers_.size(); i++)
			{
				std::vector<double> cvs_new(numnodes_, 0.0);

				if(comm_.rank() == 0)
                {
					cvs_new[mpiid_] = centers_[i];
                }

				MPI_Allreduce(MPI_IN_PLACE, &cvs_new[0], numnodes_, MPI_DOUBLE, MPI_SUM, world_);

				for(int j = 0; j < numnodes_; j++)
                {
                    worldstring_[j][i] = cvs_new[j];
                    //Represent worldstring in periodic space
                    worldstring_[j][i] = cvs[i]->GetPeriodicValue(worldstring_[j][i]); 
                }
			}
		}

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
		bool CheckEnd(const CVList& CV) 
		{
			if(maxiterator_ && iteration_ > maxiterator_)
			{
				std::cout << "System has reached max string method iterations (" << maxiterator_ << ") as specified in the input file(s)." << std::endl; 
				std::cout << "Exiting now" << std::endl; 
                PrintString(CV); //Ensure that the system prints out if it's about to exit
				world_.abort(-1);
			}

            int local_tolvalue = TolCheck();

			MPI_Allreduce(MPI_IN_PLACE, &local_tolvalue, 1, MPI_INT, MPI_LAND, world_);

			if(local_tolvalue)
			{
				std::cout << "System has converged within tolerance criteria. Exiting now" << std::endl;
                PrintString(CV); //Ensure that the system prints out if it's about to exit
				world_.abort(-1);
			}

            return true;
		}

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
		StringMethod(boost::mpi::communicator& world,
					boost::mpi::communicator& comm,
					const std::vector<double>& centers,
					unsigned int maxiterations,
					const std::vector<double> cvspring,
			 		unsigned int frequency) : 
						Method(frequency, world, comm), centers_(centers), 
						maxiterator_(maxiterations), 
						cvspring_(cvspring), iterator_(1) 
 
		{
			newcenters_.resize(centers_.size(), 0);
		}

		//! Pre-simulation hook.
		void PreSimulation(Snapshot* snapshot, const CVList& cvs) override
		{
			mpiid_ = snapshot->GetWalkerID();
			char file[1024];
			sprintf(file, "node-%04d.log",mpiid_);
		 	stringout_.open(file);

            SetSendRecvNeighbors();

			worldstring_.resize(numnodes_);
			for(auto& w : worldstring_)
				w.resize(centers_.size());
			UpdateWorldString(cvs);
			PrintString(cvs);

		}

		//! Post-integration hook.
		virtual void PostIntegration(Snapshot* snapshot, const CVList& cvs) override = 0;

		// Post-simulation hook.
		void PostSimulation(Snapshot*, const CVList&) override
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
        void SetSendRecvNeighbors()
        {
		 	std::vector<int> wiids(world_.size(), 0);

			//Set the neighbors
			recneigh_ = -1;
			sendneigh_ = -1; 

			MPI_Allgather(&mpiid_, 1, MPI_INT, &wiids[0], 1, MPI_INT, world_);
			numnodes_ = int(*std::max_element(wiids.begin(), wiids.end())) + 1;

			// Ugly for now...
			for(size_t i = 0; i < wiids.size(); i++)
			{
				if(mpiid_ == 0)
				{
					sendneigh_ = comm_.size();
					if(wiids[i] == numnodes_ - 1)
					{
						recneigh_ = i;
						break;
					}
				}
				else if (mpiid_ == numnodes_ - 1)
				{
					sendneigh_ = 0;
					if(wiids[i] == mpiid_ - 1)
					{
						recneigh_ = i;
						break;
					}
				} 
				else
				{
					if(wiids[i] == mpiid_ + 1)
					{
						sendneigh_ = i;
						break;
					}
					if(wiids[i] == mpiid_ - 1 && recneigh_ == -1)
						recneigh_ = i;
				}
			}
        }

		//! \copydoc Serializable::Serialize()
		void Serialize(Json::Value& json) const override
        {
            json["type"] = "String";

            for(size_t i = 0; i < centers_.size(); i++)
                json["centers"].append(centers_[i]);

            for(auto& t : tol_)
            	json["tolerance"].append(t);

            json["max_iterations"] = maxiterator_;

            for(auto& s : cvspring_)
            	json["ksprings"].append(s);

           json["iteration"] = iteration_; 
        }

		//! Destructor
		virtual ~StringMethod() {}
	};
}
			

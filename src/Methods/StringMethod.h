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
		std::vector<double> _centers;

		//! CV starting location values.
		std::vector<double> _newcenters;

		//! The world's strings centers for each CV.
		/*!
		 * _worldstring[node#][cv#]
		 */
		std::vector<std::vector<double> > _worldstring;

		//! The node this belongs to
		int _mpiid;

		//! Tolerance criteria for determining when to stop string (default 0 if no tolerance criteria)
		std::vector<double> _tol;

		//! Number of nodes on a string
		int _numnodes;

		//! Maximum cap on number of string method iterations performed
		unsigned int _maxiterator;

		//! Vector of spring constants.
		std::vector<double> _cvspring;

		//! The local method iterator
		unsigned int _iterator;

		//! Output stream for string data.
		std::ofstream _stringout;

		//! Neighbor to send info to.
		int _sendneigh;

		//! Neighbor to gain info from.
		int _recneigh;

		//! Store positions for starting trajectories
		std::vector<std::vector<double>> _prev_positions;

		//! Store velocities for starting trajectories
		std::vector<std::vector<double>> _prev_velocities;

		//! Store velocities for starting trajectories
		std::vector<std::vector<int>> _prev_IDs;

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
			if(_comm.rank() == 0)
			{
		        //Write node, iteration, centers of the string and current CV value to output file
		        _stringout.precision(8);
		        _stringout << _mpiid << " " << _iteration << " ";

		        for(size_t i = 0; i < _centers.size(); i++)
		            _stringout << _worldstring[_mpiid][i] << " " << CV[i]->GetValue() << " ";

		        _stringout << std::endl;
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

			if(_comm.rank() == 0)
			{
				MPI_Sendrecv(&_centers[0], _centers.size(), MPI_DOUBLE, _sendneigh, 1234,
					&(*lcv0)[0], _centers.size(), MPI_DOUBLE, _recneigh, 1234, 
					_world, &status);

				MPI_Sendrecv(&_centers[0], _centers.size(), MPI_DOUBLE, _recneigh, 4321,
			       &(*ucv0)[0], _centers.size(), MPI_DOUBLE, _sendneigh, 4321, 
			       _world, &status);
			}

			MPI_Bcast(&(*lcv0)[0],_centers.size(),MPI::DOUBLE,0,_comm);
			MPI_Bcast(&(*ucv0)[0],_centers.size(),MPI::DOUBLE,0,_comm);
		}

		//! Reparameterize the string
		/*!
		 * \param alpha_star Factor for reparametrization.
		 */
		void StringReparam(double alpha_star)
		{
			std::vector<double> alpha_star_vector(_numnodes,0.0);

			//Reparameterization
			//Alpha star is the uneven mesh, approximated as linear distance between points
			if(_comm.rank()==0)
				alpha_star_vector[_mpiid] = _mpiid == 0 ? 0 : alpha_star;

			//Gather each alpha_star into a vector 
			MPI_Allreduce(MPI::IN_PLACE, &alpha_star_vector[0], _numnodes, MPI::DOUBLE, MPI::SUM, _world);

			for(size_t i = 1; i < alpha_star_vector.size(); i++)
			    alpha_star_vector[i] += alpha_star_vector[i-1];
			
			for(size_t i = 1; i < alpha_star_vector.size(); i++)
			    alpha_star_vector[i] /= alpha_star_vector[_numnodes - 1];

			tk::spline spl; //Cubic spline interpolation

			for(size_t i = 0; i < _centers.size(); i++)
			{
				std::vector<double> cvs_new(_numnodes, 0.0);

				if(_comm.rank() == 0)
					cvs_new[_mpiid] = _centers[i];

				MPI_Allreduce(MPI::IN_PLACE, &cvs_new[0], _numnodes, MPI::DOUBLE, MPI::SUM, _world);

			    spl.set_points(alpha_star_vector, cvs_new);
			    _centers[i] = spl(_mpiid/(_numnodes - 1.0)); 
			}
		}

		//! Update the world string over MPI
		/*!
		 * \param cvs List of CVs.
		 */
		void UpdateWorldString(const CVList& cvs)
		{
			for(size_t i = 0; i < _centers.size(); i++)
			{
				std::vector<double> cvs_new(_numnodes, 0.0);

				if(_comm.rank() == 0)
                {
					cvs_new[_mpiid] = _centers[i];
                }


				MPI_Allreduce(MPI::IN_PLACE, &cvs_new[0], _numnodes, MPI::DOUBLE, MPI::SUM, _world);

				for(int j = 0; j < _numnodes; j++)
                {
                    _worldstring[j][i] = cvs_new[j];
                    //Represent worldstring in periodic space
                    _worldstring[j][i] = cvs[i]->GetPeriodicValue(_worldstring[j][i]); 
                }
			}
		}

		//! Check whether tolerance criteria has been met.
		bool TolCheck() const
		{
			if(_tol.size() == 0)
				return false;

			for(size_t i = 0; i < _centers.size(); i++)
            {
                if(fabs(_centers[i] - _worldstring[_mpiid][i]) > _tol[i])
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
			if(_maxiterator && _iteration > _maxiterator)
			{
				std::cout << "System has reached max string method iterations (" << _maxiterator << ") as specified in the input file(s)." << std::endl; 
				std::cout << "Exiting now" << std::endl; 
                PrintString(CV); //Ensure that the system prints out if it's about to exit
				_world.abort(-1);
			}

            bool local_tolvalue = TolCheck();

			MPI_Allreduce(MPI::IN_PLACE, &local_tolvalue, 1, MPI::BOOL, MPI::LAND, _world);

			if(local_tolvalue)
			{
				std::cout << "System has converged within tolerance criteria. Exiting now" << std::endl;
                PrintString(CV); //Ensure that the system prints out if it's about to exit
				_world.abort(-1);
			}

            return true;
		}

		//! Store a given snapshot.
		/*!
		 * \param snapshot Pointer to the snapshot.
		 * \param frame Index at which to store the snapshot.
		 *
		 * \note StoreSnapShot() and SetSnapshot() should be put in the snapshot
		 *       routine. Furthermore, there can be some efficiency gains made
		 *       by storing the local snapshot, and when you need to set the
		 *       snapshot, that is when you serialize.
		 */
		void StoreSnapshot(Snapshot* snapshot, int frame = 0)
		{
			auto locatoms = snapshot->GetNumAtoms();
			const auto& Pos = snapshot->GetPositions();
			const auto& Vel = snapshot->GetVelocities();
			const auto& IDs = snapshot->GetAtomIDs();

			std::vector<int> pcounts(_comm.size(), 0), mcounts(_comm.size(), 0); 
			std::vector<int> pdispls(_comm.size()+1, 0), mdispls(_comm.size()+1, 0);

			pcounts[_comm.rank()] = 3*locatoms;
			mcounts[_comm.rank()] = locatoms;

			// Reduce counts.
			MPI_Allreduce(MPI_IN_PLACE, pcounts.data(), pcounts.size(), MPI_INT, MPI_SUM, _comm);
			MPI_Allreduce(MPI_IN_PLACE, mcounts.data(), mcounts.size(), MPI_INT, MPI_SUM, _comm);

			// Compute displacements.
			std::partial_sum(pcounts.begin(), pcounts.end(), pdispls.begin() + 1);
			std::partial_sum(mcounts.begin(), mcounts.end(), mdispls.begin() + 1);

			// Re-size receiving vectors.
			_prev_positions[frame].resize(pdispls.back(), 0);
			_prev_velocities[frame].resize(pdispls.back(), 0);
			_prev_IDs[frame].resize(mdispls.back(), 0);

			std::vector<double> ptemp;
			std::vector<double> vtemp;

			for(auto& p : Pos)
			{
				ptemp.push_back(p[0]);
				ptemp.push_back(p[1]);
				ptemp.push_back(p[2]);
			}

			for(auto& v : Vel)
			{
				vtemp.push_back(v[0]);
				vtemp.push_back(v[1]);
				vtemp.push_back(v[2]);
			}

			// All-gather data.
			MPI_Allgatherv(ptemp.data(), ptemp.size(), MPI_DOUBLE, _prev_positions[frame].data(), pcounts.data(), pdispls.data(), MPI_DOUBLE, _comm);
			MPI_Allgatherv(vtemp.data(), vtemp.size(), MPI_DOUBLE, _prev_velocities[frame].data(), pcounts.data(), pdispls.data(), MPI_DOUBLE, _comm);
			MPI_Allgatherv(IDs.data(), IDs.size(), MPI_INT, _prev_IDs[frame].data(), mcounts.data(), mdispls.data(), MPI_INT, _comm);
		}

		//! Set the positions in snapshot from a previous frame.
		/*!
		 * \param snapshot Pointer to the snapshot that will be modified.
		 * \param frame Index of the frame to take the positions from.
		 */
		void SetPos(Snapshot* snapshot, int frame = 0)
		{
			auto& Pos = snapshot->GetPositions();

			for(size_t i = 0; i < _prev_IDs[frame].size(); i++)
			{
				auto localindex = snapshot->GetLocalIndex(_prev_IDs[frame][i]);
				if(localindex!= -1)
				{
					Pos[localindex][0] = _prev_positions[frame][i*3];
					Pos[localindex][1] = _prev_positions[frame][i*3 + 1];
					Pos[localindex][2] = _prev_positions[frame][i*3 + 2];
				}
			}
		}

		//! Set the velocities in a snapshot based on a previous frame.
		/*!
		 * \param snapshot Pointer to the snapshot that will be modified.
		 * \param frame Index of the frame to take the velocities from.
		 */
		void SetVel(Snapshot* snapshot, int frame = 0)
		{
			auto& Vel = snapshot->GetVelocities();

			for(size_t i = 0; i < _prev_IDs[frame].size(); i++)
			{
				auto localindex = snapshot->GetLocalIndex(_prev_IDs[frame][i]);
				if(localindex!= -1)
				{
					Vel[localindex][0] = _prev_velocities[frame][i*3];
					Vel[localindex][1] = _prev_velocities[frame][i*3 + 1];
					Vel[localindex][2] = _prev_velocities[frame][i*3 + 2];
				}
			}
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
						Method(frequency, world, comm), _centers(centers), 
						_maxiterator(maxiterations), 
						_cvspring(cvspring), _iterator(1) 
 
		{
			_newcenters.resize(_centers.size(), 0);
		}

		//! Pre-simulation hook.
		void PreSimulation(Snapshot* snapshot, const CVList& cvs) override
		{
			_mpiid = snapshot->GetWalkerID();
			char file[1024];
			sprintf(file, "node-%04d.log",_mpiid);
		 	_stringout.open(file);

            SetSendRecvNeighbors();

			_worldstring.resize(_numnodes);
			for(auto& w : _worldstring)
				w.resize(_centers.size());
			UpdateWorldString(cvs);
			PrintString(cvs);

		}

		//! Post-integration hook.
		virtual void PostIntegration(Snapshot* snapshot, const CVList& cvs) override = 0;

		// Post-simulation hook.
		void PostSimulation(Snapshot*, const CVList&) override
		{
			_stringout.close();
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
				_tol.push_back(t);
		}

		//! Communicate neighbor lists over MPI
        void SetSendRecvNeighbors()
        {
		 	std::vector<int> wiids(_world.size(), 0);

			//Set the neighbors
			_recneigh = -1;
			_sendneigh = -1; 

			MPI_Allgather(&_mpiid, 1, MPI_INT, &wiids[0], 1, MPI_INT, _world);
			_numnodes = int(*std::max_element(wiids.begin(), wiids.end())) + 1;

			// Ugly for now...
			for(size_t i = 0; i < wiids.size(); i++)
			{
				if(_mpiid == 0)
				{
					_sendneigh = _comm.size();
					if(wiids[i] == _numnodes - 1)
					{
						_recneigh = i;
						break;
					}
				}
				else if (_mpiid == _numnodes - 1)
				{
					_sendneigh = 0;
					if(wiids[i] == _mpiid - 1)
					{
						_recneigh = i;
						break;
					}
				} 
				else
				{
					if(wiids[i] == _mpiid + 1)
					{
						_sendneigh = i;
						break;
					}
					if(wiids[i] == _mpiid - 1 && _recneigh == -1)
						_recneigh = i;
				}
			}
        }

		//! \copydoc Serializable::Serialize()
		void Serialize(Json::Value& json) const override
        {
            json["type"] = "String";

            for(size_t i = 0; i < _centers.size(); i++)
                json["centers"].append(_centers[i]);

            for(auto& t : _tol)
            	json["tolerance"].append(t);

            json["max_iterations"] = _maxiterator;

            for(auto& s : _cvspring)
            	json["ksprings"].append(s);

           json["iteration"] = _iteration; 
        }

		//! Destructor
		virtual ~StringMethod() {}
	};
}
			

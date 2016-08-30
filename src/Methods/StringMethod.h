#pragma once

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
		unsigned int _mpiid;

		//! Tolerance criteria for determining when to stop string (default 0 if no tolerance criteria)
		std::vector<double> _tol;

		//! Number of nodes on a string
		unsigned int _numnodes;

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

		//! Updates the position of the string.
		virtual void StringUpdate() = 0;

		//! Helper function for calculating distances
		/*!
		 * \param x List of coordinates.
		 * \param y List of coordinates.
		 * \return Sum of distances between the x-values and y-values.
		 */
		double sqdist(const std::vector<double>& x, const std::vector<double>& y) const
		{
			double distance = 0;
			for (size_t i = 0; i < x.size(); i++)
				distance += (x[i] - y[i]) * (x[i] - y[i]);	

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
		            _stringout << _centers[i] << " " << CV[i]->GetValue() << " ";

		        _stringout << std::endl;
		    }
		}

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

		void UpdateWorldString()
		{
			for(size_t i = 0; i < _centers.size(); i++)
			{
				std::vector<double> cvs_new(_numnodes, 0.0);

				if(_comm.rank() == 0)
					cvs_new[_mpiid] = _centers[i];

				MPI_Allreduce(MPI::IN_PLACE, &cvs_new[0], _numnodes, MPI::DOUBLE, MPI::SUM, _world);

				for(size_t j = 0; j < _numnodes; j++)
					_worldstring[j][i] = cvs_new[j];
			}
		}

		//! Check whether tolerance criteria has been met.
		bool TolCheck(const CVList& cvs) const
		{
			if(_tol.size() == 0)
				return false;

			for(size_t i = 0; i < _centers.size(); i++)
				if(fabs(_centers[i] - _worldstring[_mpiid][i]) > _tol[i])
					return false;


			return true;
		}

		void CheckEnd(const CVList& cvs) const
		{
			if(_maxiterator && _iteration > _maxiterator)
			{
				std::cout << "System has reached max string method iterations (" << _maxiterator << ") as specified in the input file(s)." << std::endl;
				std::cout << "Exiting now" << std::endl;
				_world.abort(-1);
			}

			bool local_tolvalue = TolCheck(cvs);

			MPI_Allreduce(MPI::IN_PLACE, &local_tolvalue, 1, MPI::BOOL, MPI::LAND, _world);

			if(local_tolvalue)
			{
				std::cout << "System has converged within tolerance criteria. Exiting now" <<std::endl;
				_world.abort(-1);
			}
		}

	public:
		//! Constructor
		/*!
		 * \param world MPI global communicator.
		 * \param comm MPI local communicator.
		 * \param centers List of centers.
		 * \param NumNodes Number of nodes.
		 * \param maxiterator Maximum number of iterations.
		 * \param blockiterations Number of iterations per block averaging.
		 * \param tau Value of tau (default: 0.1).
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

		 	std::vector<unsigned int> wiids(_world.size(), 0);

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

			_worldstring.resize(_numnodes);
			for(auto& w : _worldstring)
				w.resize(_centers.size());
			UpdateWorldString();
			PrintString(cvs);

		}

		//! Post-integration hook.
		virtual void PostIntegration(Snapshot* snapshot, const CVList& cvs) override = 0;

		// Post-simulation hook.
		void PostSimulation(Snapshot*, const CVList&) override
		{
			_stringout.close();
		}

		void SetTolerance(std::vector<double> tol)
		{
			_tol.clear();
			for(auto& t : tol)
				_tol.push_back(t);
		}

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
			

/**
 * This file is part of
 * SSAGES - Suite for Advanced Generalized Ensemble Simulations
 *
 * Copyright 2017 Hythem Sidky <hsidky@nd.edu>
 *                Ben Sikora <bsikora906@gmail.com>
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
#include <fstream>

namespace SSAGES
{
	//! Umbrella sampling method.
	/*!
	 * Umbrella sampling method to constrain an arbitrary number of CVs at
	 * specified equilibrium distances.
	 *
	 * \ingroup Methods
	 */
	class Umbrella : public Method, public Buildable<Umbrella>
	{
	private:
		//! Vector of spring constants.
		std::vector<double> kspring_;

		//! Vector of equilibrium distances.
		std::vector<double> centers0_, centers1_;

		//! Amount of time over which to scale centers.
		int time_;

		//! File name
		std::string filename_;

		//! Log every n time steps
		int logevery_;

		//! Output stream for umbrella data.
		std::ofstream umbrella_;

		double GetCurrentCenter(int iteration, unsigned i)
		{
			// We are at the end.
			if(iteration >= time_) return centers1_[i];

			// Scale linearly.
			return (centers1_[i] - centers0_[i])/time_*iteration + centers0_[i]; 
		}

		//! Print umbrella values.
		/*!
		 * \param cvs List of CVs.
		 * \param iteration Current iteration.
		 */
		void PrintUmbrella(const CVList& cvs, uint iteration);

	public:
		//! Constructor.
		/*!
		 * \param world MPI global communicator.
		 * \param comm MPI local communicator.
		 * \param kspring List of spring constants.
		 * \param centers List of spring equilibrium positions.
		 * \param name Filename.
		 * \param frequency Frequency with which this method is applied.
		 *
		 * Create instance of umbrella with spring constants "kspring", and
		 * centers "centers". Note the sizes of the vectors should be
		 * commensurate with the number of CVs.
		 */
		Umbrella(const MPI_Comm& world,
				 const MPI_Comm& comm,
				 const std::vector<double>& kspring,
				 const std::vector<double>& centers,
				 std::string name,
				 unsigned int frequency) : 
		Method(frequency, world, comm), kspring_(kspring), centers0_(centers),
		centers1_(centers), time_(0), filename_(name), logevery_(1)
		{
		}

		//! Constructor.
		/*!
		 * \param world MPI global communicator.
		 * \param comm MPI local communicator.
		 * \param kspring List of spring constants.
		 * \param centers0 List of starting spring equilibrium positions.
		 * \param centers1 List of ending spring equilibrium positions.
		 * \param timesteps Number of timesteps over which to scale centers.
		 * \param name Filename.
		 * \param frequency Frequency with which this method is applied.
		 *
		 * Create instance of umbrella with spring constants "kspring", and
		 * centers "centers". Note the sizes of the vectors should be
		 * commensurate with the number of CVs.
		 */
		Umbrella(const MPI_Comm& world,
				 const MPI_Comm& comm,
				 const std::vector<double>& kspring,
				 const std::vector<double>& centers0,
				 const std::vector<double>& centers1,
				 int timesteps,
				 std::string name,
				 unsigned int frequency) : 
		Method(frequency, world, comm), kspring_(kspring), centers0_(centers0),
		centers1_(centers1), time_(timesteps), filename_(name), logevery_(1)
		{
		}

		//! Pre-simulation hook.
		/*!
		 * \param snapshot Simulation snapshot.
		 * \param cvs List of CVs.
		 */
		void PreSimulation(Snapshot* snapshot, const CVList& cvs) override;

		//! Post-integration hook.
		/*!
		 * \param snapshot Simulation snapshot.
		 * \param cvs List of CVs.
		 */
		void PostIntegration(Snapshot* snapshot, const CVList& cvs) override;

		//! Post-simulation hook.
		/*!
		 * \param snapshot Simulation snapshot.
		 * \param cvs List of CVs.
		 */
		void PostSimulation(Snapshot* snapshot, const CVList& cvs) override;

		//! Set how often to log
		/*!
		 * \param iter New value for logging interval.
		 */
		void SetLogStep(const int iter)
		{
			logevery_ = iter;
		}

		//! \copydoc Buildable::Build()
		static Umbrella* Build(const Json::Value& json, 
		                       const MPI_Comm& world,
		                       const MPI_Comm& comm,
					           const std::string& path);

		//! \copydoc Serializable::Serialize()
		/*!
		 * \warning The serialization is not implemented yet.
		 */
		void Serialize(Json::Value& json) const override
		{
			json["type"] = "Umbrella";
			for(auto& k : kspring_)
				json["ksprings"].append(k);

			if(time_ != 0 )
			{
				for(auto& c : centers0_)
					json["centers0"].append(c);
				
				for(auto& c : centers1_)
					json["centers1"].append(c);

				json["timesteps"] = time_;
			}
			else
			{			
				for(auto& c : centers0_)
					json["centers"].append(c);
			}

			json["file name"] = filename_;
			json["log every"] = logevery_;
		}

	};
}

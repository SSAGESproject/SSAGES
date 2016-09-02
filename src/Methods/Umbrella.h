/**
 * This file is part of
 * SSAGES - Suite for Advanced Generalized Ensemble Simulations
 *
 * Copyright 2016 Hythem Sidky <hsidky@nd.edu>
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
#include "../CVs/CollectiveVariable.h"
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
	class Umbrella : public Method
	{
	private:
		//! Vector of spring constants.
		std::vector<double> _kspring;

		//! Vector of equilibrium distances.
		std::vector<double> _centers;

		//! File name
		std::string _filename;

		//! Log every n time steps
		int _logevery;

		//! Output stream for umbrella data.
		std::ofstream _umbrella;

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
		Umbrella(boost::mpi::communicator& world,
				 boost::mpi::communicator& comm,
				 const std::vector<double>& kspring,
				 const std::vector<double>& centers,
				 std::string name,
				 unsigned int frequency) : 
		Method(frequency, world, comm), _kspring(kspring), _centers(centers),
		_filename(name), _logevery(1)
		{}

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
		
		//! Print umbrella values.
		/*!
		 * \param cvs List of CVs.
		 */
		void PrintUmbrella(const CVList& cvs);

		//! Set how often to log
		/*!
		 * \param iter New value for logging interval.
		 */
		void SetLogStep(const int iter)
		{
			_logevery = iter;
		}

		//! \copydoc Serializable::Serialize()
		/*!
		 * \warning The serialization is not implemented yet.
		 */
		void Serialize(Json::Value& json) const override
		{
			json["type"] = "Umbrella";
			for(auto& k : _kspring)
				json["ksprings"].append(k);

			for(auto& c : _centers)
				json["centers"].append(c);

			json["file name"] = _filename;

			json["iteration"] = _iteration;

			json["log every"] = _logevery;

		}

	};
}

/**
 * This file is part of
 * SSAGES - Suite for Advanced Generalized Ensemble Simulations
 *
 * Copyright 2016 Hythem Sidky <hsidky@nd.edu>
 *                Ben Sikora <bsikora906@gmail.com>
 *                Emre Sevgen <sesevgen@uchicago.edu>
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

#include "CollectiveVariable.h"
#include "Validator/ObjectRequirement.h"
#include "Drivers/DriverException.h"
#include "Snapshot.h"
#include "schema.h"
#include <array>
#include <cmath>

namespace SSAGES
{
	//! Collective variable on the torsion angles.
	/*!
	 * Collective variable on an proper dihedral. This will return the angle
	 * between two planes as defined in \cite VANSCHAIK1993751. Singularities
	 * are avoided as described in \cite BLONDEL19961132.
	 *
	 * \ingroup CVs
	 */
	class TorsionalCV : public CollectiveVariable
	{
	private:
		//! Vector of 4 atom ID's of interest.
		Label atomids_;

		//! If \c True, use periodic boundary conditions.
		bool periodic_;

	public:
		//! Constructor.
		/*!
		 * \param atomid1 ID of the first atom defining the dihedral angle.
		 * \param atomid2 ID of the second atom defining the dihedral angle.
		 * \param atomid3 ID of the third atom defining the dihedral angle.
		 * \param atomid4 ID of the forth atom defining the dihedral angle.
		 * \param periodic If \c True consider periodic boundary conditions.
		 *
		 * Construct an dihedral CV.
		 *
		 * \todo Bounds needs to be an input and periodic boundary conditions
		 */
		TorsionalCV(int atomid1, int atomid2, int atomid3, int atomid4, bool periodic) : 
		atomids_({atomid1, atomid2, atomid3, atomid4}), periodic_(periodic)
		{
		}

		//! Initialize necessary variables.
		/*!
		 * \param snapshot Current simulation snapshot.
		 */
		void Initialize(const Snapshot& snapshot) override
		{
			using std::to_string;
			
			std::vector<int> found;
			snapshot.GetLocalIndices(atomids_, &found);
			int nfound = found.size();
			MPI_Allreduce(MPI_IN_PLACE, &nfound, 1, MPI_INT, MPI_SUM, snapshot.GetCommunicator());
			
			if(nfound != 4)
				throw BuildException({
					"TorsionalCV: Expected to find " + 
					to_string(4) + 
					" atoms, but only found " + 
					to_string(nfound) + "."
				});		
		}

		//! Evaluate the CV.
		/*!
		 * \param snapshot Current simulation snapshot.
		 */
		void Evaluate(const Snapshot& snapshot) override
		{
			// Get data from snapshot. 
			auto n = snapshot.GetNumAtoms();
			const auto& pos = snapshot.GetPositions();
			auto& comm = snapshot.GetCommunicator();

			// Initialize gradient.
			std::fill(grad_.begin(), grad_.end(), Vector3{0,0,0});
			grad_.resize(n, Vector3{0,0,0});

			Vector3 xi{0, 0, 0}, xj{0, 0, 0}, xk{0, 0, 0}, xl{0, 0, 0};

			auto iindex = snapshot.GetLocalIndex(atomids_[0]);
			auto jindex = snapshot.GetLocalIndex(atomids_[1]);
			auto kindex = snapshot.GetLocalIndex(atomids_[2]);
			auto lindex = snapshot.GetLocalIndex(atomids_[3]);

			if(iindex != -1) xi = pos[iindex];
			if(jindex != -1) xj = pos[jindex];
			if(kindex != -1) xk = pos[kindex];
			if(lindex != -1) xl = pos[lindex];

			// By performing a reduce, we actually collect all. This can 
			// be converted to a more intelligent allgater on rank then bcast. 
			MPI_Allreduce(MPI_IN_PLACE, xi.data(), 3, MPI_DOUBLE, MPI_SUM, comm);
			MPI_Allreduce(MPI_IN_PLACE, xj.data(), 3, MPI_DOUBLE, MPI_SUM, comm);
			MPI_Allreduce(MPI_IN_PLACE, xk.data(), 3, MPI_DOUBLE, MPI_SUM, comm);
			MPI_Allreduce(MPI_IN_PLACE, xl.data(), 3, MPI_DOUBLE, MPI_SUM, comm);
			
			//Calculate pertinent vectors
			auto F = snapshot.ApplyMinimumImage(xi - xj);
			auto G = snapshot.ApplyMinimumImage(xj - xk);
			auto H = snapshot.ApplyMinimumImage(xl - xk);
			auto A = F.cross(G);
			auto B = H.cross(G);
			
			auto y = B.cross(A).dot(G.normalized());
			auto x = A.dot(B);

			val_ = atan2(y, x);

			auto Zed = F.dot(G.normalized())/A.dot(A); 
			auto Ned = H.dot(G.normalized())/B.dot(B);

			Vector3 gradi{0,0,0}, gradl{0,0,0}; 
			if(iindex != -1) gradi = -G.norm()*A/A.dot(A);
			if(lindex != -1) gradl = G.norm()*B/B.dot(B);
			MPI_Allreduce(MPI_IN_PLACE, gradi.data(), 3, MPI_DOUBLE, MPI_SUM, comm);
			MPI_Allreduce(MPI_IN_PLACE, gradl.data(), 3, MPI_DOUBLE, MPI_SUM, comm);
			if(iindex != -1) grad_[iindex] = gradi;
			if(lindex != -1) grad_[lindex] = gradl;
			if(jindex != -1) grad_[jindex] = Zed*A - Ned*B - gradi;
			if(kindex != -1) grad_[kindex] = Ned*B  - Zed*A - gradl;
		}

		//! Return value taking periodic boundary conditions into account
		/*!
		 * \param Location Input value.
		 * \return Wrapped or unwrapped input value depending on whether
		 *         periodic boundaries are used.
		 *
		 * If periodic boundaries are used, this function wraps the input
		 * value into the range (-pi, pi). Otherwise the input value is
		 * returned unmodified.
		 */
		double GetPeriodicValue(double Location) const override
		{
			if(!periodic_)
				return Location;

			int n = (int)(Location/(2.0*M_PI));
			double PeriodicLocation;

			PeriodicLocation = Location - n*M_PI;
			if(PeriodicLocation < -M_PI)
				PeriodicLocation += 2.0*M_PI;
			else if (PeriodicLocation > M_PI)
				PeriodicLocation -= 2.0*M_PI;

			return PeriodicLocation;
		}

		//! Get difference taking periodic boundary conditions into account.
		/*!
		 * \param Location Input value
		 * \return Wrapped or unwrapped difference depending on whether periodic
		 *         boundaries are used.
		 *
		 * If periodic boundaries are used, this function calculates the
		 * difference and wraps the result into the range (-pi, pi). Otherwise,
		 * the simple difference is returned.
		 */
		double GetDifference(const double Location) const override
		{
			double PeriodicDiff = val_ - Location;

			if(!periodic_)
				return PeriodicDiff;

			PeriodicDiff = GetPeriodicValue(PeriodicDiff);

			if(PeriodicDiff > M_PI)
				PeriodicDiff -= 2.0*M_PI;
			else if(PeriodicDiff < -M_PI)
				PeriodicDiff += 2.0*M_PI;

			return PeriodicDiff;
		}

        //! Returns the minimum image of a CV based on the input location.
        /*!
		 * \param Location Value against which the minimum image is calculated.
		 * \return Minimum image of the CV 
		 *
         * Takes the input location and applies the periodic boundary conditions to return a minimum image
         * of the CV.
		 */
		double GetMinimumImage(const double Location) const override
		{
			double PeriodicDiff = val_ - Location;

			if(!periodic_)
				return val_;	

			if(PeriodicDiff > M_PI)
				return (val_ - 2.0*M_PI);
			else if(PeriodicDiff < -M_PI)
				return (val_ + 2.0*M_PI);	

			return val_;
		}

		//! \copydoc CollectiveVariable::BuildCV()
		static TorsionalCV* Build(const Json::Value& json, const std::string& path)
		{
			Json::ObjectRequirement validator;
			Json::Value schema;
			Json::Reader reader;

			reader.parse(JsonSchema::TorsionalCV, schema);
			validator.Parse(schema, path);

			// Validate inputs.
			validator.Validate(json, path);
			if(validator.HasErrors())
				throw BuildException(validator.GetErrors());

			std::vector<int> atomids;
			for(auto& s : json["atom_ids"])
				atomids.push_back(s.asInt());

			auto periodic = json.get("periodic", true).asBool();
			return new TorsionalCV(atomids[0], atomids[1], atomids[2], atomids[3], periodic);
		}
	};
}

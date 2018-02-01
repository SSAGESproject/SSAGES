/**
 * This file is part of
 * SSAGES - Suite for Advanced Generalized Ensemble Simulations
 *
 * Copyright 2017 Michael Webb <xmwebb@gmail.com>
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

namespace SSAGES
{
	//! Collective variable is a Rouse mode for a polymer chain comprised of N particle groups
	/*!
	 *  This CV returns the value for the p'th Rouse mode, computed as
	 *  Xp(t) ~ (1/N) sum_{i=1}^{N} ri(t)*cos[pi(n-0.5)p/N],
	 *  where N is the number of particle groups, p is the mode index, ri is the center-of-mass position of
	 *  a collection of atoms comprising the i'th bead in the N-bead polymer chain
	 *  \ingroup CVs
	 */
	class RouseModeCV : public CollectiveVariable
	{
	private:
		std::vector<Label>     groups_; //!< vector of groups of indices to define the particle groups
		std::vector<double>     massg_; //!< vector of the total mass for each particle group
		size_t                      p_; //!< index of mode of interest as CV
		size_t                      N_; //!< number of Rouse beads in polymer chain
		Vector3                    xp_; //!< 3d vector for containing vectorial rouse amplitude
		std::vector<Vector3>        r_; //!< vector of coordinate positions for each bead
	
	public:
		//! Basic Constructor for Rouse Mode CV
		/*!
		 * \param groups - vector of vector of atom IDs, each group comprising a bead in the Rouse chain
		 * \param p      - index for relevant Rouse mode
		 */
		RouseModeCV(const std::vector<Label>& groups, int p) : 
		groups_(groups), p_(p), N_( groups_.size())
		{ massg_.resize( groups_.size(),0.0 ); }

		//! Helper function to determine masses of each group
		/*! Note: this here assumes that the masses of each group are not changing during
 		 *        the simulation, which is likely typical...
 		 * \param snapshot Current simulation snapshot.
 		 * \param groups - vector of vector of atom IDs, each group comprising a bead in the Rouse chain
 		 */ 
		void setMasses(const std::vector<Label>& groups, const Snapshot& snapshot ) {
			// Compute mass of each group and store to massg_
			Label listi;
			for (size_t i = 0; i < N_; ++i) {
				listi     = groups[i];     // MW: can this be used directly in Total Mass
				Label idi;
				snapshot.GetLocalIndices(listi, &idi);
				massg_[i] = snapshot.TotalMass(idi);
			} 
		}

		//! Initialize necessary variables.
		/*!
		 * \param snapshot Current simulation snapshot.
		 */
		void Initialize(const Snapshot& snapshot) override
		{
			using std::to_string;

			// Check for valid p_
			if (p_ > groups_.size())
				throw std::invalid_argument(
					"RouseModeCV: Expected to find p to be less than " + 
					to_string(groups_.size()) +" but found p = " + 
					to_string(p_) 
				);			
			// Check for valid groups
			for (size_t j = 0; j < groups_.size(); ++j) {
				auto nj = groups_[j].size();
				// Make sure atom IDs in the group are somewhere
				std::vector<int> found(nj,0);
				for (size_t i = 0; i < nj; ++i) {
					if(snapshot.GetLocalIndex(groups_[j][i]) != -1)
						found[i] = 1;
				}

				MPI_Allreduce(MPI_IN_PLACE, found.data(), nj, MPI_INT, MPI_SUM, snapshot.GetCommunicator());
				unsigned njtot = std::accumulate(found.begin(), found.end(), 0, std::plus<int>());
				
				if(njtot != nj)
					throw std::invalid_argument(
						"RouseModeCV: Expected to find " + 
						to_string(nj) + 
						" atoms in group " + to_string(j) +", but only found " + 
						to_string(njtot) + "."
					);			

			}
			// Set the masses of each particle group in massg_
			this->setMasses(groups_, snapshot);
		}

		//! Evaluate the CV.
		/*!
		 * \param snapshot Current simulation snapshot.
		 */
		void Evaluate(const Snapshot& snapshot) override
		{

			// Get data from snapshot.
			auto          ntot = snapshot.GetNumAtoms(); // total number of atoms
			const auto& masses = snapshot.GetMasses();   // mass of each atom

			// Initialize working variables
			double ppi_n = p_ * M_PI / N_;  // constant
			xp_.fill(0.0);  	// vectorial Rouse mode
			r_.resize(N_);		// position vector for beads in Rouse chain (unwrapped)
			grad_.resize(ntot, Vector3{0,0,0});  // gradient set to 0.0 for all atoms
			std::fill(grad_.begin(), grad_.end(), Vector3{0,0,0});

			// Iterate over all N_ atom groups and compute the center of mass for each
			std::vector<Vector3>     rcom; // vector of COM positions 
			for (size_t i = 0; i < N_; ++i) {
				Label idi;	// list of indices
				snapshot.GetLocalIndices(groups_[i], &idi);
				rcom.push_back(snapshot.CenterOfMass(idi,massg_[i])); // center of mass for group
			}

			// Now compute differences vectors between the neighboring beads
			// accumulate displacements to reconstruct unwrapped polymer chain
			// for simplicity, we consider the first bead to be the reference position
			// in all snapshots
			r_[0] = rcom[0];
			for (size_t i = 1; i < N_; ++i) {
				Vector3 dri = snapshot.ApplyMinimumImage(rcom[i] - rcom[i-1]);
				r_[i] = r_[i-1] + dri;  // r_i = r_{i-1} + (r_i - r_{i-1})
			}

			// Determine the value of the Rouse coordinate
			// Xp(t) = 1/sqrt(N) * sum_{i=1}^{N} ri, p = 0
			// Xp(t) = sqrt(2/N) * sum_{i=1}^{N} ri * cos[p*pi/N*(i-0.5)], p = 1,...,N-1
			// Note: this solution is valid for homogeneous friction
			xp_.fill(0.0);
			for (size_t i = 0; i < N_; ++i) {
				xp_ += r_[i]*cos ( ppi_n * (i+0.5) );
			}
			xp_ /= sqrt(N_);
			if ( p_ != 0 ) xp_ *= sqrt(2.0) ;

			// Compute Rouse vector norm as the CV
			// CV = sqrt(Xp*Xp), Xp = (Xp1,Xp2,Xp3)
			val_ = xp_.norm();

			// Now perform gradient operation
			// dCV/dxjd = (Xpd/CV)*(c/N)**0.5*sum_i=1^N cos[p*pi(i-0.5)/N] mj/Mi *delta_j({i})
			// delta_j({i}) = 1, if j in {i}, 0 otherwise
			Vector3   gradpcon  = xp_ / sqrt(N_) / val_;
			if ( p_ != 0) gradpcon *= sqrt(2.0);         // (Xpd/CV)*(c/N)**0.5
			for (size_t i = 0; i < N_; ++i) {
				Label idi;	// list of indices
				snapshot.GetLocalIndices(groups_[i], &idi);
				// go over each atom in the group and add to its gradient
				// Note: performance tradeoff here. All gradient elements have a common factor of
				// Xpd/CV*sqrt(c/N) = prefactor, with c = 2 if p != 0
				// this could be factored out, but if ntot >> number of atoms in groups
				// then it won't be worth it to post multiply all gradient terms... 
				// Could also go over loop again after to do the multiplication, but 
				// that is troublesome if atom ids appear in multiple groups for some reason
				double cosval = cos(ppi_n*(i+0.5)) / massg_[i]; // cos[p*pi(i-0.5)/N] / Mi
				for (auto& id : idi) {
					grad_[id] +=  gradpcon* cosval * masses[id];
				}
			}

		}

		//! \copydoc CollectiveVariable::BuildCV()
		static RouseModeCV* Build(const Json::Value& json, const std::string& path)
		{
			Json::ObjectRequirement validator;
			Json::Value schema;
			Json::Reader reader;

			reader.parse(JsonSchema::RouseModeCV, schema);
			validator.Parse(schema, path);

			// Validate inputs.
			validator.Validate(json, path);
			if(validator.HasErrors())
				throw BuildException(validator.GetErrors());
			
			std::vector<Label> groups;
			for (auto& group : json["groups"]) 
			{
				groups.push_back({});
				for(auto& id : group) 
					groups.back().push_back(id.asInt());
			}

			auto mode = json.get("mode",0).asInt();
			return new RouseModeCV( groups, mode);	
		}
	};
}

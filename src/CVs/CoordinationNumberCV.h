/**
 * This file is part of
 * SSAGES - Suite for Advanced Generalized Ensemble Simulations
 *
 * Copyright 2017 Hythem Sidky <hsidky@nd.edu>
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

 #include "CVs/CollectiveVariable.h"
 #include "Drivers/DriverException.h"
 
 namespace SSAGES
 {
	class SwitchingFunction : public Serializable
	{
	private: 
		double d0_; //!< Minimum linear shift value. 
		double r0_; //!< Cutoff distance. 
		int n_, m_; //!< Exponents of the switching function which controll the stiffness. 
	public:
		SwitchingFunction(double d0, double r0, int n, int m) : 
		d0_(d0), r0_(r0), n_(n), m_(m) {}

		//! Evaluate the switching function.
		/*!
		 * \param rij distance between two atoms. 
		 * \param df Reference to variable which will store the gradient.
		 * 
		 * \return value of switching function. 
		 */
		double Evaluate(double rij, double& df) const
		{
			const auto xarg = (rij - d0_)/r0_;
			const auto xn = std::pow(xarg, n_);
			const auto xm = std::pow(xarg, m_);
			const auto f = (1-xn)/(1-xm);
			
			df = f/(d0_-rij)*(n_*xn/(1-xn)+m_*xm/(xm-1));
			return 	f;
		}

		//! Build SwitchingFunction from JSON value. 
		/*!
		 * \param json JSON value node. 
		 * 
		 * \return Pointer to new SwitchingFunction.
		 */
		static SwitchingFunction* Build(Json::Value& json)
		{
			return new SwitchingFunction(
			                json["d0"].asDouble(), 
							json["r0"].asDouble(), 
							json["n"].asInt(), 
							json["m"].asInt()
						);
		}

		//! Serialize this CV for restart purposes.
		/*!
		 * \param json JSON value
		 */
		void Serialize(Json::Value& json) const override
		{
			json["d0"] = d0_;
			json["r0"] = r0_;
			json["n"] = n_;
			json["m"] = m_;
		}
	};

    //! Collective variable on coordination number between two groups of atoms.
	/*!
	 * Collective variable on coordination number between two groups of atoms. 
     * To ensure a continuously differentiable function, there are various 
     * switching functions from which to choose from. 
	 *
	 * \ingroup CVs
	 */
	class CoordinationNumberCV : public CollectiveVariable
	{
	private:
		Label group1_; //!< IDs of the first group of atoms. 
		Label group2_; //!< IDs of the second group of atoms. 
		SwitchingFunction sf_; //!< Switching function used for coordination number.

	public:
		//! Constructor.
		/*!
		 * \param group1 IDs of the first group of atoms. 
		 * \param group2 IDs of the second group of atoms. 
		 *
		 * Construct a CoordinationNumberCV.
		 *
		 */    
		CoordinationNumberCV(const Label& group1, const Label& group2, SwitchingFunction&& sf) : 
		group1_(group1), group2_(group2), sf_(sf)
		{
		}

		//! Initialize necessary variables.
		/*!
		 * \param snapshot Current simulation snapshot.
		 */
		void Initialize(const Snapshot& snapshot) override
		{
			using std::to_string;

			auto n1 = group1_.size(), n2 = group2_.size();

			// Make sure atom ID's are on at least one processor. 
			std::vector<int> found1(n1, 0), found2(n2, 0);
			for(size_t i = 0; i < n1; ++i)
			{
				if(snapshot.GetLocalIndex(group1_[i]) != -1)
					found1[i] = 1;
			}
			
			for(size_t i = 0; i < n2; ++i)
			{
				if(snapshot.GetLocalIndex(group2_[i]) != -1)
					found2[i] = 1;
			}

			MPI_Allreduce(MPI_IN_PLACE, found1.data(), n1, MPI_INT, MPI_SUM, snapshot.GetCommunicator());
			MPI_Allreduce(MPI_IN_PLACE, found2.data(), n2, MPI_INT, MPI_SUM, snapshot.GetCommunicator());
	
			unsigned ntot1 = std::accumulate(found1.begin(), found1.end(), 0, std::plus<int>());
			unsigned ntot2 = std::accumulate(found2.begin(), found2.end(), 0, std::plus<int>());
			if(ntot1 != n1)
				throw BuildException({
					"CoordinationNumberCV: Expected to find " + 
					to_string(n1) + 
					" atoms in group 1, but only found " + 
					to_string(ntot1) + "."
				});			

			if(ntot2 != n2)
				throw BuildException({
					"CoordinationNumberCV: Expected to find " + 
					to_string(n2) + 
					" atoms in group 2, but only found " + 
					to_string(ntot2) + "."
				});			
		}
		
		//! Evaluate the CV.
		/*!
		 * \param snapshot Current simulation snapshot.
		 */
		void Evaluate(const Snapshot& snapshot) override
		{
			// Get local atom indices.
			std::vector<int> idx1, idx2;
			snapshot.GetLocalIndices(group1_, &idx1);
			snapshot.GetLocalIndices(group2_, &idx2);

			// Get data from snapshot. 
			auto n = snapshot.GetNumAtoms();
			auto& atomids = snapshot.GetAtomIDs();
			auto& positions = snapshot.GetPositions();

			// Initialize gradient.
			std::fill(grad_.begin(), grad_.end(), Vector3{0,0,0});
			grad_.resize(n, Vector3{0,0,0});

			// The nastiness begins. We essentially need to compute 
			// pairwise distances between the atoms. For now, let's 
			// gather the atomic coordinates of group2_
			// on all relevant processors. Can be improved later. 
			auto& comm = snapshot.GetCommunicator();
			std::vector<int> pcounts(comm.size(), 0), icounts(comm.size(), 0); 
			std::vector<int> pdispls(comm.size()+1, 0), idispls(comm.size()+1, 0); 
			pcounts[comm.rank()] = 3*idx2.size();
			icounts[comm.rank()] = idx2.size();

			// Reduce counts.
			MPI_Allreduce(MPI_IN_PLACE, pcounts.data(), pcounts.size(), MPI_INT, MPI_SUM, comm);
			MPI_Allreduce(MPI_IN_PLACE, icounts.data(), icounts.size(), MPI_INT, MPI_SUM, comm);

			std::partial_sum(pcounts.begin(), pcounts.end(), pdispls.begin() + 1);
			std::partial_sum(icounts.begin(), icounts.end(), idispls.begin() + 1);
			
			// Fill up mass and position vectors.
			std::vector<double> pos(3*idx2.size(), 0);
			std::vector<int> ids(idx2.size(), 0);
			for(size_t i = 0; i < idx2.size(); ++i)
			{
				auto& idx = idx2[i];
				auto& p = positions[idx];
				pos[3*i+0] = p[0];
				pos[3*i+1] = p[1];
				pos[3*i+2] = p[2];
				ids[i] = atomids[idx];
			}

			// Create receiving vectors.
			std::vector<double> gpos(pdispls.back(), 0);
			std::vector<int> gids(idispls.back(), 0);
			
			// Gather data.
			MPI_Allgatherv(pos.data(), pos.size(), MPI_DOUBLE, gpos.data(), pcounts.data(), pdispls.data(), MPI_DOUBLE, comm);
			MPI_Allgatherv(ids.data(), ids.size(), MPI_INT, gids.data(), icounts.data(), idispls.data(), MPI_INT, comm);

			val_ = 0;
			double df = 0.;
			for(auto& i : idx1)
			{
				auto& p = positions[i];
				for(size_t j = 0; j < group2_.size(); ++j)
				{
					Vector3 rij = {p[0] - gpos[3*j], p[1] - gpos[3*j+1], p[2] - gpos[3*j+2]};
					snapshot.ApplyMinimumImage(rij);
					auto r = rij.norm();
					val_ +=  sf_.Evaluate(r, df);

					grad_[i] += df*rij/r;

					auto lidx = snapshot.GetLocalIndex(gids[j]);
					if(lidx != -1)
						grad_[lidx] -= df*rij/r;
				}
			}
			
			// Reduce val. 
			MPI_Allreduce(MPI_IN_PLACE, &val_, 1, MPI_DOUBLE, MPI_SUM, comm);
		}
		
		//! Serialize this CV for restart purposes.
		/*!
		 * \param json JSON value
		 */
		void Serialize(Json::Value& json) const override
		{
			json["type"] = "CoordinationNumber";
		}
	};
 }
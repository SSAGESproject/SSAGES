/**
 * This file is part of
 * SSAGES - Suite for Advanced Generalized Ensemble Simulations
 *
 * Copyright 2017 Michael Quevillon <mquevill@nd.edu>
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
#include "Validator/ObjectRequirement.h"
#include "Drivers/DriverException.h"
#include "Snapshot.h"
#include "schema.h"
 
namespace SSAGES
{
	class GaussianFunction : public Serializable
	{
	private: 
		double mu_; //!< Center of Gaussian.
		double sigma_; //!< Width of Gaussian 
	public:
		GaussianFunction(double mu, double sigma) : 
		mu_(mu), sigma_(sigma) {}

		//! Evaluate the Gaussian function.
		/*!
			* \param rij distance between two atoms. 
			* \param df Reference to variable which will store the gradient.
			* 
			* \return value of Gaussian function. 
			*/
		double Evaluate(double rij, double& df) const
		{
			const auto dx = (rij - mu_)/sigma_;
			const auto f = exp( - dx*dx/2.);
			const auto pre = - dx/sigma_;
			
			df = pre * f;
			return 	f;
		}

		//! Build GaussianFunction from JSON value. 
		/*!
			* \param json JSON value node. 
			* 
			* \return Pointer to new GaussianFunction.
			*/
		static GaussianFunction Build(const Json::Value& json)
		{
			return GaussianFunction(
						json["mu"].asDouble(), 
						json["sigma"].asDouble()
					);
		}

		//! Serialize this CV for restart purposes.
		/*!
			* \param json JSON value
			*/
		void Serialize(Json::Value& json) const override
		{
			json["mu"] = mu_;
			json["sigma"] = sigma_;
		}
	};

	//! Collective variable on sum of Gaussian kernel applied to pairwise distances.
	/*!
		* Collective variable on sum of Gaussian kernel applied to pairwise distances.
		* If Gaussians parameters are chosen judiciously, can emulate counting nearest
		* neighbors. (For example, a steep Gaussian centered at contact distance should
		* return a number close to the number of nearest neighbors.) However, this CV
		* is a continuous function, so it may not return the exact integer value.
		*
		* \ingroup CVs
		*/
	class NearestNeighborsCV : public CollectiveVariable, public Buildable<NearestNeighborsCV>
	{
	private:
		Label group1_; //!< IDs of the (only) group of atoms. 
		GaussianFunction gf_; //!< Gaussian function used for calculating

	public:
		//! Constructor.
		/*!
			* \param group1 IDs of the (only) group of atoms.
			*
			* Construct a NearestNeighborsCV.
			*
			*/    
		NearestNeighborsCV(const Label& group1, GaussianFunction&& gf) : 
		group1_(group1), gf_(gf)
		{
		}

		//! Initialize necessary variables.
		/*!
			* \param snapshot Current simulation snapshot.
			*/
		void Initialize(const Snapshot& snapshot) override
		{
			using std::to_string;

			auto n1 = group1_.size();

			// Make sure atom ID's are on at least one processor. 
			std::vector<int> found1(n1, 0);
			for(size_t i = 0; i < n1; ++i)
			{
				if(snapshot.GetLocalIndex(group1_[i]) != -1)
					found1[i] = 1;
			}

			MPI_Allreduce(MPI_IN_PLACE, found1.data(), n1, MPI_INT, MPI_SUM, snapshot.GetCommunicator());

			unsigned ntot1 = std::accumulate(found1.begin(), found1.end(), 0, std::plus<int>());
			if(ntot1 != n1)
				throw BuildException({
					"NearestNeighborsCV: Expected to find " + 
					to_string(n1) + 
					" atoms in group 1, but only found " + 
					to_string(ntot1) + "."
				});
		}
		
		//! Evaluate the CV.
		/*!
			* \param snapshot Current simulation snapshot.
			*/
		void Evaluate(const Snapshot& snapshot) override
		{
			// Get local atom indices.
			std::vector<int> idx1;
			snapshot.GetLocalIndices(group1_, &idx1);

			// Get data from snapshot. 
			auto n = snapshot.GetNumAtoms();
			auto& atomids = snapshot.GetAtomIDs();
			auto& positions = snapshot.GetPositions();

			// Initialize gradient.
			std::fill(grad_.begin(), grad_.end(), Vector3{0,0,0});
			grad_.resize(n, Vector3{0,0,0});
			boxgrad_ = Matrix3::Zero();

			// Need to compute pairwise distances between all the atoms.
			auto& comm = snapshot.GetCommunicator();
			std::vector<int> pcounts(comm.size(), 0), icounts(comm.size(), 0); 
			std::vector<int> pdispls(comm.size()+1, 0), idispls(comm.size()+1, 0); 
			pcounts[comm.rank()] = 3*idx1.size();
			icounts[comm.rank()] = idx1.size();

			// Reduce counts.
			MPI_Allreduce(MPI_IN_PLACE, pcounts.data(), pcounts.size(), MPI_INT, MPI_SUM, comm);
			MPI_Allreduce(MPI_IN_PLACE, icounts.data(), icounts.size(), MPI_INT, MPI_SUM, comm);

			std::partial_sum(pcounts.begin(), pcounts.end(), pdispls.begin() + 1);
			std::partial_sum(icounts.begin(), icounts.end(), idispls.begin() + 1);
			
			// Fill up mass and position vectors.
			std::vector<double> pos(3*idx1.size(), 0);
			std::vector<int> ids(idx1.size(), 0);
			for(size_t i = 0; i < idx1.size(); ++i)
			{
				auto& idx = idx1[i];
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
				for(auto& j : idx1)
				{
					if(i != j)
					{
						auto& q = positions[j];
						Vector3 rij = {p[0] - q[0], p[1] - q[1], p[2] - q[2]};
						snapshot.ApplyMinimumImage(&rij);
						auto r = rij.norm();
						val_ +=  gf_.Evaluate(r, df);

						grad_[i] += df*rij/r;
						boxgrad_ += df*rij/r*rij.transpose();
					}
				}
				grad_[i] /= 2.0; // Correct for double counting
			}
			val_ /= 2.0; // Correct for double counting

			// Reduce val. 
			MPI_Allreduce(MPI_IN_PLACE, &val_, 1, MPI_DOUBLE, MPI_SUM, comm);
		}

		static NearestNeighborsCV* Construct(const Json::Value& json, const std::string& path)
		{
			Json::ObjectRequirement validator;
			Json::Value schema;
			Json::Reader reader;

			reader.parse(JsonSchema::NearestNeighborsCV, schema);
			validator.Parse(schema, path);

			// Validate inputs.
			validator.Validate(json, path);
			if(validator.HasErrors())
				throw BuildException(validator.GetErrors());
			
			std::vector<int> group1;
			
			for(auto& s : json["group1"])
				group1.push_back(s.asInt());
			
			return new NearestNeighborsCV(group1, GaussianFunction::Build(json["gaussian"]));
		}

		//! Serialize this CV for restart purposes.
		/*!
			* \param json JSON value
			*/
		void Serialize(Json::Value& json) const override
		{
			json["type"] = "NearestNeighbors";
		}
	};
}

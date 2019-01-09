/**
 * This file is part of
 * SSAGES - Software Suite for Advanced General Ensemble Simulations
 *
 * Copyright 2017 Hythem Sidky <hsidky@nd.edu>
 *                Michael Quevillon <mquevill@nd.edu>
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
#include "Utility/PairwiseKernel.h"
#include "Snapshot.h"
#include "schema.h"
 
namespace SSAGES
{
	//! Generalized collective variable based on pairwise properties of atoms.
	/*!
	 * Collective variable on pairwise properties between two groups of atoms. 
	 * To ensure generality of usage, there are various pairwise kernel functions
	 * from which to choose. 
	 *
	 * \ingroup CVs
	 */
	class PairwiseCV : public CollectiveVariable
	{
	private:
		Label group1_; //!< IDs of the first group of atoms. 
		Label group2_; //!< IDs of the second group of atoms. 
		PairwiseKernel* pk_; //!< Pairwise kernel function used for CV.

	public:
		//! Constructor.
		/*!
		 * \param group1 IDs of the first group of atoms. 
		 * \param group2 IDs of the second group of atoms. 
		 * \param pk pairwise kernel (function) to use.
		 *
		 * Construct a PairwiseCV.
		 *
		 */    
		PairwiseCV(const Label& group1, const Label& group2, PairwiseKernel* pk) : 
		group1_(group1), group2_(group2), pk_(pk)
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
					"Pairwise: Expected to find " + 
					to_string(n1) + 
					" atoms in group 1, but only found " + 
					to_string(ntot1) + "."
				});			

			if(ntot2 != n2)
				throw BuildException({
					"Pairwise: Expected to find " + 
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
			auto& comm = snapshot.GetCommunicator();

			// Initialize gradient.
			std::fill(grad_.begin(), grad_.end(), Vector3{0,0,0});
			grad_.resize(n, Vector3{0,0,0});
			boxgrad_ = Matrix3::Zero();

			// Fill eigen matrix with coordinate data. 
			std::vector<double> pos1(3*idx1.size()), pos2(3*idx2.size());
			std::vector<int> id1(idx1.size()), id2(idx2.size()); 
			for(size_t i = 0; i < idx1.size(); ++i)
			{
				pos1[3*i+0] = positions[idx1[i]][0]; 
				pos1[3*i+1] = positions[idx1[i]][1]; 
				pos1[3*i+2] = positions[idx1[i]][2];

				id1[i] = atomids[idx1[i]];
			}

			for(size_t i = 0; i < idx2.size(); ++i)
			{
				pos2[3*i+0] = positions[idx2[i]][0]; 
				pos2[3*i+1] = positions[idx2[i]][1]; 
				pos2[3*i+2] = positions[idx2[i]][2];

				id2[i] = atomids[idx2[i]];
			}

			// Gather across all procs. 
			pos1 = mxx::allgatherv(pos1.data(), pos1.size(), comm);
			pos2 = mxx::allgatherv(pos2.data(), pos2.size(), comm);
			id1 = mxx::allgatherv(id1.data(), id1.size(), comm);
			id2 = mxx::allgatherv(id2.data(), id2.size(), comm);

			// Create Eigen map for ease of use.
			using Map = Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor>>;
			Map gpos1(pos1.data(), group1_.size(), 3), gpos2(pos2.data(), group2_.size(), 3);

			val_ = 0;
			double df = 0.;
			// Compute gradient and value.
			for(size_t i = 0; i < group1_.size(); ++i)
			{
				auto t1 = id1[i];
				const auto& pi = gpos1.row(i);
				for(size_t j = 0; j < group2_.size(); ++j)
				{
					auto t2 = id2[j];
					const auto& pj = gpos2.row(j);

					// Skip identical pairs. 
					if(t1 == t2)
						continue;
					
					// Compute distance and pairwise function. 
					Vector3 rij = pi - pj;
					snapshot.ApplyMinimumImage(&rij);
					auto r = rij.norm();
					val_ +=  pk_->Evaluate(r, df);

					// Get local indices and sum gradients.
					auto lid1 = snapshot.GetLocalIndex(t1);
					if(lid1 != -1)
						grad_[lid1] += df*rij/r;
					
					auto lid2 = snapshot.GetLocalIndex(t2);
					if(lid2 != -1)
						grad_[lid2] -= df*rij/r;

					// Only sum box gradient on a single proc.
					if(comm.rank() == 0)
						boxgrad_ += df*rij/r*rij.transpose();
				}
			}
		}

		//! \copydoc CollectiveVariable::BuildCV()
		static PairwiseCV* Build(const Json::Value& json, const std::string& path)
		{
			Json::ObjectRequirement validator;
			Json::Value schema;
			Json::CharReaderBuilder rbuilder;
			Json::CharReader* reader = rbuilder.newCharReader();

			reader->parse(JsonSchema::PairwiseCV.c_str(),
			              JsonSchema::PairwiseCV.c_str() + JsonSchema::PairwiseCV.size(),
			              &schema, nullptr);
			validator.Parse(schema, path);

			// Validate inputs.
			validator.Validate(json, path);
			if(validator.HasErrors())
				throw BuildException(validator.GetErrors());
			
			std::vector<int> group1, group2;
			
			for(auto& s : json["group1"])
				group1.push_back(s.asInt());

			for(auto& s : json["group2"])
				group2.push_back(s.asInt());
			
			return new PairwiseCV(group1, group2, PairwiseKernel::Build(json["kernel"], path));
		}

		//! Destructor
		~PairwiseCV()
		{
			delete pk_;
		}
	};
}

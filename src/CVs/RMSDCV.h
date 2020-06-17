/**
 * This file is part of
 * SSAGES - Software Suite for Advanced General Ensemble Simulations
 *
 * Copyright 2016 Yamil Colon <yamilcolon2015@u.northwestern.edu>
 *                Hythem Sidky <hsidky@nd.edu>
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
#include "Utility/ReadFile.h"
#include <array>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

namespace SSAGES
{
	//! Collective variable to calculate root mean square displacement.
	/*!
	 * RMSD calculation performed as reported in
	 * Coutsias, E. A., Seok, C., and Dill, K. A., "Using Quaternions to Calculate RMSD", J. Comput. Chem. 25: 1849-1857, 2004
	 */

	class RMSDCV : public CollectiveVariable
	{
	private:
		//! IDs of the atoms used for RMSD calculation
		Label atomids_;

		//! Name of model structure.
		std::string xyzfile_;

		//! Store reference structure coordinates.
		std::vector<Vector3> refcoord_;

		//! Center of mass of reference.
		Vector3 COMref_;

	public:
		//! Constructor.
		/*!
		 * \param atomids IDs of the atoms defining RMSD.
		 * \param xyzfile String determining the reference file.
		 * \param use_range If \c True Use range of atoms defined by the two atoms in atomids.
		 *
		 * Construct a RMSD CV.
		 *
		 * \todo Bounds needs to be an input and periodic boundary conditions
		 */
		RMSDCV(const Label& atomids, std::string xyzfile, bool use_range = false) :
		atomids_(atomids), xyzfile_(xyzfile)
		{
			if(use_range)
			{
				if(atomids.size() != 2)
				{
					throw BuildException({
						"RMSDCV: With use_range, expected 2 atoms, but found " +
						std::to_string(atomids.size()) + "."
					});
				}

				atomids_.clear();
				if(atomids[0] >= atomids[1])
				{
					throw BuildException({
						"RMSDCV: Atom range must be strictly increasing."
					});
				}
				for(int i = atomids[0]; i <= atomids[1]; ++i)
				{
					atomids_.push_back(i);
				}
			}
		}

		//! Initialize necessary variables.
		/*!
		 * \param snapshot Current simulation snapshot.
		 */
		void Initialize(const Snapshot& snapshot) override
		{
			const auto& masses = snapshot.GetMasses();
			auto n = atomids_.size();

			// Make sure atom IDs are on at least one processor.
			std::vector<int> found(n, 0);
			for(size_t i = 0; i < n; ++i)
			{
				if(snapshot.GetLocalIndex(atomids_[i]) != -1) found[i] = 1;
			}

			MPI_Allreduce(MPI_IN_PLACE, found.data(), n, MPI_INT, MPI_SUM, snapshot.GetCommunicator());
			unsigned ntot = std::accumulate(found.begin(), found.end(), 0, std::plus<int>());
			if(ntot != n)
			{
				throw BuildException({
					"RMSDCV: Expected to find " + std::to_string(n) +
					" atoms, but only found " + std::to_string(ntot) + "."
				});
			}

			Label idx;
			snapshot.GetLocalIndices(atomids_, &idx);

			std::vector<std::array<double,4>> xyzinfo = ReadFile::ReadXYZ(xyzfile_);

			refcoord_.resize(snapshot.GetNumAtoms(), Vector3{0,0,0});

			// Loop through atom positions
			for(size_t i = 0; i < xyzinfo.size(); ++i)
			{
				int id = snapshot.GetLocalIndex(xyzinfo[i][0]);
				if(id == -1) continue; // Atom not found
				for(size_t j = 0; j < atomids_.size(); ++j)
				{
					if(atomids_[j] == xyzinfo[i][0])
					{
						refcoord_[id][0] = xyzinfo[i][1];
						refcoord_[id][1] = xyzinfo[i][2];
						refcoord_[id][2] = xyzinfo[i][3];
					}
				}
			}

			Vector3 mass_pos_prod_ref = Vector3::Zero();
			double total_mass = 0;

			for(auto& i : idx)
			{
				mass_pos_prod_ref += masses[i]*refcoord_[i];
				total_mass += masses[i];
			}
			MPI_Allreduce(MPI_IN_PLACE, mass_pos_prod_ref.data(), mass_pos_prod_ref.size(), MPI_DOUBLE, MPI_SUM, snapshot.GetCommunicator());
			MPI_Allreduce(MPI_IN_PLACE, &total_mass, 1, MPI_DOUBLE, MPI_SUM, snapshot.GetCommunicator());

			COMref_ = mass_pos_prod_ref/total_mass;
		}

		//! Evaluate the CV.
		/*!
		 * \param snapshot Current simulation snapshot.
		 */
		void Evaluate(const Snapshot& snapshot) override
		{
			// Get local atom indices.
			Label idx;
			snapshot.GetLocalIndices(atomids_, &idx);

			// Get data from snapshot.
			const auto& masses = snapshot.GetMasses();
			const auto& pos = snapshot.GetPositions();
			double masstot = snapshot.TotalMass(idx);

			// Initialize gradient.
			std::fill(grad_.begin(), grad_.end(), Vector3{0,0,0});
			grad_.resize(snapshot.GetNumAtoms(), Vector3{0,0,0});

			// Calculate center of mass
			Vector3 COM_ = snapshot.CenterOfMass(idx);

			// Build correlation matrix
			Matrix3 R = Matrix3::Zero();

			Vector3 diff, diff_ref;
			double part_RMSD = 0;

			for(auto& i : idx)
			{
				diff = snapshot.ApplyMinimumImage(pos[i] - COM_);
				diff_ref = refcoord_[i] - COMref_; // Reference has no "box" or minimum image

				R.noalias() += masses[i]*diff*diff_ref.transpose();

				part_RMSD += masses[i]*(diff.squaredNorm() + diff_ref.squaredNorm());
			}
			MPI_Allreduce(MPI_IN_PLACE, R.data(), R.size(), MPI_DOUBLE, MPI_SUM, snapshot.GetCommunicator());
			MPI_Allreduce(MPI_IN_PLACE, &part_RMSD, 1, MPI_DOUBLE, MPI_SUM, snapshot.GetCommunicator());
			R /= masstot;
			part_RMSD /= masstot;

			Eigen::Matrix4d F;
			F(0,0) = R(0,0) + R(1,1) + R(2,2);
			F(1,0) = R(1,2) - R(2,1);
			F(2,0) = R(2,0) - R(0,2);
			F(3,0) = R(0,1) - R(1,0);

			F(0,1) = F(1,0);
			F(1,1) = R(0,0) - R(1,1) - R(2,2);
			F(2,1) = R(0,1) + R(1,0);
			F(3,1) = R(0,2) + R(2,0);

			F(0,2) = F(2,0);
			F(1,2) = F(2,1);
			F(2,2) = -R(0,0) + R(1,1) - R(2,2);
			F(3,2) = R(1,2) + R(2,1);

			F(0,3) = F(3,0);
			F(1,3) = F(3,1);
			F(2,3) = F(3,2);
			F(3,3) = -R(0,0) - R(1,1) + R(2,2);

			//Find eigenvalues
			Eigen::SelfAdjointEigenSolver<Eigen::Matrix4d> es(F);
			auto max_lambda = es.eigenvalues().real()[3];

			// Calculate RMSD
			val_ = sqrt(part_RMSD - (2*max_lambda));

			// If RMSD is zero, we have nothing to do.
			if(val_ == 0) return;

			// Calculate gradient
			auto ev = es.eigenvectors().real().col(3);
			auto q = Eigen::Quaterniond(ev(0),ev(1),ev(2),ev(3));
			Matrix3 RotMatrix = q.toRotationMatrix();

			for(auto& i : idx)
			{
				auto rotref = RotMatrix.transpose()*(refcoord_[i] - COMref_);
				grad_[i] = masses[i]/masstot*(1-masses[i]/masstot)*(snapshot.ApplyMinimumImage(pos[i] - COM_) - rotref)/val_;
			}
		}

		//! \copydoc CollectiveVariable::BuildCV()
		static RMSDCV* Build(const Json::Value& json, const std::string& path)
		{
			Json::ObjectRequirement validator;
			Json::Value schema;
			Json::CharReaderBuilder rbuilder;
			Json::CharReader* reader = rbuilder.newCharReader();

			reader->parse(JsonSchema::RMSDCV.c_str(),
			              JsonSchema::RMSDCV.c_str() + JsonSchema::RMSDCV.size(),
			              &schema, NULL);
			validator.Parse(schema, path);

			//Validate inputs
			validator.Validate(json, path);
			if(validator.HasErrors())
					throw BuildException(validator.GetErrors());

			std::vector<int> atom_ids;
			for(auto& s : json["atom_ids"])
				atom_ids.push_back(s.asInt());
			auto reference = json["reference"].asString();
			auto use_range = json.get("use_range",false).asBool();

			return new RMSDCV(atom_ids, reference, use_range);
		}
	};
}

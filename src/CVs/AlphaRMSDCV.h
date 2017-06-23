/**
 * This file is part of
 * SSAGES - Suite for Advanced Generalized Ensemble Simulations
 *
 * Copyright 2017 Ashley Guo <azguo@uchicago.edu>
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
#include "Utility/ReadBackbone.h"

namespace SSAGES
{
	//! Collective variable to measure alpha helix secondary structure
	/*!
	 * Following treatment in Pietrucci and Laio, "A Collective Variable for
	 * the Efficient Exploration of Protein Beta-Sheet Structures: Application
	 * to SH3 and GB1", JCTC, 2009, 5(9): 2197-2201.
	 *
	 * Check blocks of six consecutive protein residues for RMSD from
	 * reference "ideal" alpha helix structure.
	 */

	class AlphaRMSDCV : public CollectiveVariable
	{
	private:

		//!< Residue IDs for secondary structure calculation
		std::vector<int> resids_;

		//!< Atom IDs for secondary structure calculation: backbone of resids_
		std::vector<int> atomids_;

		//!< Name of pdb reference for system
		std::string refpdb_;

		//!< Coordinates for reference structure
		std::vector<Vector3> refalpha_;

		//!< Length unit conversion: convert 1 nm to your internal MD units (ex. if using angstroms use 10)
		double unitconv_;

	public:
		//! Constructor.
		/*!
		 * \param resids IDs of residues for calculating secondary structure
		 * \param refpdb String of pdb filename with atom and residue indices.
		 *
		 * Construct an AlphaRMSD CV -- description.
		 *
		 * \todo Bounds needs to be an input and periodic boundary conditions ?
		 */
		AlphaRMSDCV(std::vector<int> resids, std::string refpdb, double unitconv) :
		resids_(resids), refpdb_(refpdb), unitconv_(unitconv)
		{
			if(resids_.size() != 2 ){
				std::cout << "AlphaRMSDCV: Input must designate range of residues with 2 residue numbers." << std::endl;
				exit(0);
			}

			resids_.clear();

			if(resids[0] >= resids[1]){
				std::cout << "AlphaRMSDCV: Input must list lower residue index first: please reverse residue range." << std::endl;
				exit(0);
			} else if(resids[1] - resids[0] < 5) {
				std::cout << "AlphaRMSDCV: Residue range must span at least 6 residues for alpha helix calculation." << std::endl;
				exit(0);
			}

			std::cout << "AlphaRMSDCV: Calculating alpha helix character from residue " << resids[0] << " to " << resids[1]  << "." << std::endl;

			for(unsigned int i = resids[0]; i <= resids[1]; i++){
				resids_.push_back(i);
			}
		}

		// Initialize variables
		void Initialize(const Snapshot& snapshot) override
		{
			atomids_ = ReadBackbone::GetPdbBackbone(refpdb_, resids_);

			// reference 'ideal' alpha helix structure, in nanometers
			refalpha_.push_back(unitconv_ * Vector3{ .0733,  .0519,  .5298 }); // N
			refalpha_.push_back(unitconv_ * Vector3{ .1763,  .0810,  .4301 }); // CA
			refalpha_.push_back(unitconv_ * Vector3{ .3166,  .0543,  .4881 }); // CB
			refalpha_.push_back(unitconv_ * Vector3{ .1527, -.0045,  .3053 }); // C
			refalpha_.push_back(unitconv_ * Vector3{ .1646,  .0436,  .1928 }); // O
			refalpha_.push_back(unitconv_ * Vector3{ .1180, -.1312,  .3254 }); // N
			refalpha_.push_back(unitconv_ * Vector3{ .0924, -.2203,  .2126 }); // CA
			refalpha_.push_back(unitconv_ * Vector3{ .0650, -.3626,  .2626 }); // CB
			refalpha_.push_back(unitconv_ * Vector3{-.0239, -.1711,  .1261 }); // C
			refalpha_.push_back(unitconv_ * Vector3{-.0190, -.1815,  .0032 }); // O
			refalpha_.push_back(unitconv_ * Vector3{-.1280, -.1172,  .1891 }); // N
			refalpha_.push_back(unitconv_ * Vector3{-.2416, -.0661,  .1127 }); // CA
			refalpha_.push_back(unitconv_ * Vector3{-.3548, -.0217,  .2056 }); // CB
			refalpha_.push_back(unitconv_ * Vector3{-.1964,  .0529,  .0276 }); // C
			refalpha_.push_back(unitconv_ * Vector3{-.2364,  .0659, -.0880 }); // O
			refalpha_.push_back(unitconv_ * Vector3{-.1130,  .1391,  .0856 }); // N
			refalpha_.push_back(unitconv_ * Vector3{-.0620,  .2565,  .0148 }); // CA
			refalpha_.push_back(unitconv_ * Vector3{ .0228,  .3439,  .1077 }); // CB
			refalpha_.push_back(unitconv_ * Vector3{ .0231,  .2129, -.1032 }); // C
			refalpha_.push_back(unitconv_ * Vector3{ .0179,  .2733, -.2099 }); // O
			refalpha_.push_back(unitconv_ * Vector3{ .1028,  .1084, -.0833 }); // N
			refalpha_.push_back(unitconv_ * Vector3{ .1872,  .0593, -.1919 }); // CA
			refalpha_.push_back(unitconv_ * Vector3{ .2850, -.0462, -.1397 }); // CB
			refalpha_.push_back(unitconv_ * Vector3{ .1020,  .0020, -.3049 }); // C
			refalpha_.push_back(unitconv_ * Vector3{ .1317,  .0227, -.4224 }); // O
			refalpha_.push_back(unitconv_ * Vector3{-.0051, -.0684, -.2696 }); // N
			refalpha_.push_back(unitconv_ * Vector3{-.0927, -.1261, -.3713 }); // CA
			refalpha_.push_back(unitconv_ * Vector3{-.1933, -.2219, -.3074 }); // CB
			refalpha_.push_back(unitconv_ * Vector3{-.1663, -.0171, -.4475 }); // C
			refalpha_.push_back(unitconv_ * Vector3{-.1916, -.0296, -.5673 }); // O
		}

		// Evaluate the CV
		void Evaluate(const Snapshot& snapshot) override
		{
			// need atom positions for all atoms in atomids_
			const auto& pos = snapshot.GetPositions();
			//const auto& image_flags = snapshot.GetImageFlags();
			std::vector<int> groupidx;
			snapshot.GetLocalIndices(atomids_, &groupidx);	// get correct local atom indices

			double rmsd, dist_norm, dxgrouprmsd;
			Vector3 dist_xyz, dist_ref;
			std::vector<Vector3> refxyz;
			std::vector< std::vector< Vector3 > > deriv(30, std::vector<Vector3>(30, Vector3{0,0,0}));

			std::fill(grad_.begin(), grad_.end(), Vector3{0,0,0});
			grad_.resize(snapshot.GetNumAtoms(), Vector3{0,0,0});
			val_ = 0.0;

			unsigned int resgroups = resids_.size() - 5;

			// for each set of 6 residues
			for(size_t i = 0; i < resgroups; i++){

				// clear temp rmsd calculation
				rmsd = 0.0;

				// load refxyz with the correct 30 reference atoms
				std::fill(refxyz.begin(), refxyz.end(), Vector3{0,0,0});
				refxyz.resize(30, Vector3{0,0,0});
				for(size_t j = 0; j < 30; j++){
					refxyz[j] = pos[groupidx[5 * i + j]];
				}

				// sum over 435 pairs in refxyz and refalpha_ to calculate CV
				for(size_t j = 0; j < 29; j++){
					for(size_t k = j + 1; k < 30; k++){
						dist_xyz = refxyz[j] - refxyz[k];
						dist_ref = refalpha_[j] - refalpha_[k]; // could be precalculated
						dist_norm = dist_xyz.norm() - dist_ref.norm();
						rmsd += dist_norm * dist_norm;
						deriv[j][k] = dist_xyz * dist_norm / dist_xyz.norm();
					}
				}

				rmsd = pow(rmsd / 435, 0.5) / 0.1;
				val_ += (1 - pow(rmsd, 8.0)) / (1 - pow(rmsd, 12.0));

				// is ther ea better WAY
				dxgrouprmsd = 8E8 * (5E3 * pow(rmsd, 11.0) + pow(rmsd, 7.0));
				dxgrouprmsd /= 1E8 * pow(rmsd, 8.0) + 1E4 * pow(rmsd, 4.0) + 1;
				dxgrouprmsd /= 1E8 * pow(rmsd, 8.0) + 1E4 * pow(rmsd, 4.0) + 1;
				dxgrouprmsd *= -2.0 / 435;  // derivative has prefactor of -2 / (30 choose 2)

				for(size_t j = 0; j < 29; j++){
					for(size_t k = j + 1; k < 30; k++){
						grad_[groupidx[5 * i + j]] += dxgrouprmsd * deriv[j][k];
						grad_[groupidx[5 * i + k]] -= dxgrouprmsd * deriv[j][k];
					}
				}
			}
		}

		//! Serialize this CV for restart purposes.
		/*!
		 * \param json JSON value
		 */
		virtual void Serialize(Json::Value& json) const override
		{
			json["type"] = "AlphaRMSD";
			json["reference"] = refpdb_;
			for(size_t i=0; i < resids_.size(); ++i)
				json["residue_ids"].append(resids_[i]);
			json["length_unit"] = unitconv_;
		}
	};
}

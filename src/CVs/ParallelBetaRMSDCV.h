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
	//! Collective variable to measure parallel beta sheet secondary structure
	/*!
	 * Following treatment in Pietrucci and Laio, "A Collective Variable for
	 * the Efficient Exploration of Protein Beta-Sheet Structures: Application
	 * to SH3 and GB1", JCTC, 2009, 5(9): 2197-2201.
	 *
	 * Check 2 blocks of 3 consecutive protein residues for RMSD from
	 * reference "ideal" parallel beta sheet structure.
	 */

	class ParallelBetaRMSDCV : public CollectiveVariable
	{
	private:

		//!< Residue IDs for secondary structure calculation
		std::vector<int> resids_;

		//!< Atom IDs for secondary structure calculation: backbone of resids_
		//std::vector<int> atomids_;
		std::vector< std::vector<std::string> > atomids_;

		//!< Name of pdb reference for system
		std::string refpdb_;

		//!< Coordinates for reference structure
		std::vector<Vector3> refalpha_;

		//!< Length unit conversion: convert 1 nm to your internal MD units (ex. if using angstroms use 10)
		double unitconv_;

		//!< mode
		int mode_;

	public:
		//! Constructor.
		/*!
		 * \param resids IDs of residues for calculating secondary structure
		 * \param refpdb String of pdb filename with atom and residue indices.
		 *
		 * \todo Bounds needs to be an input and periodic boundary conditions ?
		 */
		ParallelBetaRMSDCV(std::vector<int> resids, std::string refpdb, double unitconv, int mode) :
		resids_(resids), refpdb_(refpdb), unitconv_(unitconv), mode_(mode)
		{
			if(resids_.size() != 2 ){
				std::cout << "ParallelBetaRMSDCV: Input must designate range of residues with 2 residue numbers." << std::endl;
				exit(0);
			}

			resids_.clear();

			if(resids[0] >= resids[1]){
				std::cout << "ParallelBetaRMSDCV: Input must list lower residue index first: please reverse residue range." << std::endl;
				exit(0);
			} else if(resids[1] - resids[0] < 5) {
				std::cout << "ParallelBetaRMSDCV: Residue range must span at least 6 residues for secondary structure calculation." << std::endl;
				exit(0);
			}

			std::cout << "ParallelBetaRMSDCV: Calculating parallel beta sheet character from residue " << resids[0] << " to " << resids[1]  << "." << std::endl;

			for(unsigned int i = resids[0]; i <= resids[1]; i++){
				resids_.push_back(i);
			}
		}

		// Initialize variables
		void Initialize(const Snapshot& snapshot) override
		{
			atomids_ = ReadBackbone::GetPdbBackbone(refpdb_, resids_);

			// reference structure for parallel beta sheet, values in nanometers
			refalpha_.push_back(unitconv_ * Vector3{-.1439, -.5122, -.1144}); // N    residue i
			refalpha_.push_back(unitconv_ * Vector3{-.0816, -.3803, -.1013}); // CA
			refalpha_.push_back(unitconv_ * Vector3{ .0099, -.3509, -.2206}); // CB
			refalpha_.push_back(unitconv_ * Vector3{-.1928, -.2770, -.0952}); // C
			refalpha_.push_back(unitconv_ * Vector3{-.2991, -.2970, -.1551}); // O
			refalpha_.push_back(unitconv_ * Vector3{-.1698, -.1687, -.0215}); // N    residue i+1
			refalpha_.push_back(unitconv_ * Vector3{-.2681, -.0613, -.0143}); // CA
			refalpha_.push_back(unitconv_ * Vector3{-.3323, -.0477,  .1267}); // CB
			refalpha_.push_back(unitconv_ * Vector3{-.1984,  .0681, -.0574}); // C
			refalpha_.push_back(unitconv_ * Vector3{-.0807,  .0921, -.0273}); // O
			refalpha_.push_back(unitconv_ * Vector3{-.2716,  .1492, -.1329}); // N    residue i+2
			refalpha_.push_back(unitconv_ * Vector3{-.2196,  .2731, -.1883}); // CA
			refalpha_.push_back(unitconv_ * Vector3{-.2263,  .2692, -.3418}); // CB
			refalpha_.push_back(unitconv_ * Vector3{-.2989,  .3949, -.1433}); // C
			refalpha_.push_back(unitconv_ * Vector3{-.4214,  .3989, -.1583}); // O
			refalpha_.push_back(unitconv_ * Vector3{ .2464, -.4352,  .2149}); // N    residue h
			refalpha_.push_back(unitconv_ * Vector3{ .3078, -.3170,  .1541}); // CA
			refalpha_.push_back(unitconv_ * Vector3{ .3398, -.3415,  .0060}); // CB
			refalpha_.push_back(unitconv_ * Vector3{ .2080, -.2021,  .1639}); // C
			refalpha_.push_back(unitconv_ * Vector3{ .0938, -.2178,  .1225}); // O
			refalpha_.push_back(unitconv_ * Vector3{ .2525, -.0886,  .2183}); // N    residue h+1
			refalpha_.push_back(unitconv_ * Vector3{ .1692,  .0303,  .2346}); // CA
			refalpha_.push_back(unitconv_ * Vector3{ .1541,  .0665,  .3842}); // CB
			refalpha_.push_back(unitconv_ * Vector3{ .2420,  .1410,  .1608}); // C
			refalpha_.push_back(unitconv_ * Vector3{ .3567,  .1733,  .1937}); // O
			refalpha_.push_back(unitconv_ * Vector3{ .1758,  .1976,  .0600}); // N    residue h+2
			refalpha_.push_back(unitconv_ * Vector3{ .2373,  .2987, -.0238}); // CA
			refalpha_.push_back(unitconv_ * Vector3{ .2367,  .2527, -.1720}); // CB
			refalpha_.push_back(unitconv_ * Vector3{ .1684,  .4331, -.0148}); // C
			refalpha_.push_back(unitconv_ * Vector3{ .0486,  .4430, -.0415}); // O
		}

		// Evaluate the CV
		void Evaluate(const Snapshot& snapshot) override
		{
			// need atom positions for all atoms in atomids_
			const auto& pos = snapshot.GetPositions();
			std::vector<int> groupidx;
			std::vector<int> double_atomids(atomids_[0].size());
			std::transform(atomids_[0].begin(), atomids_[0].end(), double_atomids.begin(), [](std::string val) {return std::stod(val);});
			snapshot.GetLocalIndices(double_atomids, &groupidx);

			unsigned int resgroups = resids_.size() - 2;
			double rmsd, dist_norm, dxgrouprmsd;
			Vector3 dist_xyz, dist_ref;
			std::vector<Vector3> refxyz;
			std::vector< std::vector< Vector3 > > deriv(30, std::vector<Vector3>(30, Vector3{0,0,0}));

			std::fill(grad_.begin(), grad_.end(), Vector3{0,0,0});
			grad_.resize(snapshot.GetNumAtoms(), Vector3{0,0,0});
			val_ = 0.0;

			for(size_t i = 0; i < resgroups - 3; i++){
				for(size_t j = i + 3; j < resgroups; j++){
					// mode: 0 for all, 1 for inter, 2 for intra
					if((mode_ == 0) || (mode_ == 1 && atomids_[1][5 * j] != atomids_[1][5 * i]) || (mode_ == 2 && atomids_[1][5 * j] == atomids_[1][5 * i])){
						rmsd = 0.0;
						std::fill(refxyz.begin(), refxyz.end(), Vector3{0,0,0});
						refxyz.resize(30, Vector3{0,0,0});
						for(size_t k = 0; k < 15; k++){
							refxyz[k] = pos[groupidx[5 * i + k]];
						}
						for(size_t k = 0; k < 15; k++){
							refxyz[k + 15] = pos[groupidx[5 * j + k]];
						}

						// sum over all pairs to calculate CV
						for(size_t k = 0; k < 29; k++){
							for(size_t h = k + 1; h < 30; h++){
								dist_xyz = refxyz[k] - refxyz[h];
								dist_ref = refalpha_[k] - refalpha_[h];
								dist_norm = dist_xyz.norm() - dist_ref.norm();
								rmsd += dist_norm * dist_norm;
								deriv[k][h] = dist_xyz * dist_norm / dist_xyz.norm();
							}
						}

						rmsd = pow(rmsd / 435, 0.5) / 0.1;
						val_ += (1 - pow(rmsd, 8.0)) / (1 - pow(rmsd, 12.0));

						dxgrouprmsd = pow(rmsd, 11.0) + pow(rmsd, 7.0);
						dxgrouprmsd /= pow(rmsd, 8.0) + pow(rmsd, 4.0) + 1;
						dxgrouprmsd /= pow(rmsd, 8.0) + pow(rmsd, 4.0) + 1;
						dxgrouprmsd *= -40. / 435;

						for(size_t k = 0; k < 15; k++){
							for(size_t h = k + 1; h < 15; h++){
								grad_[groupidx[5 * i + k]] += dxgrouprmsd * deriv[k][h];
								grad_[groupidx[5 * i + h]] -= dxgrouprmsd * deriv[k][h];
							}
							for(size_t h = 0; h < 15; h++){
								grad_[groupidx[5 * i + k]] += dxgrouprmsd * deriv[k][h+15];
								grad_[groupidx[5 * j + h]] -= dxgrouprmsd * deriv[k][h+15];
							}
						}
						for(size_t k = 0; k < 14; k++){
							for(size_t h = k + 1; h < 15; h++){
								grad_[groupidx[5 * j + k]] += dxgrouprmsd * deriv[k+15][h+15];
								grad_[groupidx[5 * j + h]] -= dxgrouprmsd * deriv[k+15][h+15];
							}
						}
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
			json["type"] = "ParallelBetaRMSD";
			json["reference"] = refpdb_;
			for(size_t i=0; i < resids_.size(); ++i)
				json["residue_ids"].append(resids_[i]);
			json["length_unit"] = unitconv_;
		}
	};
}

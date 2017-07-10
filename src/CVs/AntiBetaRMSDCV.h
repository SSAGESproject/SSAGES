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
#include "Validator/ObjectRequirement.h"
#include "Drivers/DriverException.h"
#include "Snapshot.h"
#include "schema.h"
#include "Utility/ReadBackbone.h"
#include <stdexcept>

namespace SSAGES
{
	//! Collective variable to measure antiparallel beta secondary structure
	/*!
	 * Following treatment in Pietrucci and Laio, "A Collective Variable for
	 * the Efficient Exploration of Protein Beta-Sheet Structures: Application
	 * to SH3 and GB1", JCTC, 2009, 5(9): 2197-2201.
	 *
	 * Check 2 blocks of 3 consecutive protein residues for RMSD from
	 * reference "ideal" antiparallel beta sheet structure.
	 */

	class AntiBetaRMSDCV : public CollectiveVariable, public Buildable<AntiBetaRMSDCV>
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
		AntiBetaRMSDCV(std::vector<int> resids, std::string refpdb, double unitconv, int mode) :
		resids_(resids), refpdb_(refpdb), unitconv_(unitconv), mode_(mode)
		{
			if(resids_.size() != 2 ){
				throw std::invalid_argument("AntiBetaRMSDCV: Input must designate range of residues with 2 residue numbers.");
			}

			resids_.clear();

			if(resids[0] >= resids[1]){
				throw std::invalid_argument("AntiBetaRMSDCV: Input must list lower residue index first: please reverse residue range.");
			} else if(resids[1] - resids[0] < 5){
				throw std::invalid_argument("AntiBetaRMSDCV: Residue range must span at least 6 residues for secondary structure calculation.");
			}

			std::cout << "AntiBetaRMSDCV: Calculating antiparallel beta sheet character from residue " << resids[0] << " to " << resids[1]  << "." << std::endl;

			for(unsigned int i = resids[0]; i <= resids[1]; i++){
				resids_.push_back(i);
			}
		}

		// Initialize variables
		void Initialize(const Snapshot& snapshot) override
		{
			atomids_ = ReadBackbone::GetPdbBackbone(refpdb_, resids_);

			// reference 'ideal' antiparallel beta sheet structure in nanometers
			// Reference structure taken from values used in Plumed 2.2.0
			refalpha_.push_back(unitconv_ * Vector3{ .2263, -.3795,  .1722}); // N    i
			refalpha_.push_back(unitconv_ * Vector3{ .2493, -.2426,  .2263}); // CA
			refalpha_.push_back(unitconv_ * Vector3{ .3847, -.1838,  .1761}); // CB
			refalpha_.push_back(unitconv_ * Vector3{ .1301, -.1517,  .1921}); // C
			refalpha_.push_back(unitconv_ * Vector3{ .0852, -.1504,  .0739}); // O
			refalpha_.push_back(unitconv_ * Vector3{ .0818, -.0738,  .2917}); // N    i+1
			refalpha_.push_back(unitconv_ * Vector3{-.0299,  .0243,  .2748}); // CA
			refalpha_.push_back(unitconv_ * Vector3{-.1421, -.0076,  .3757}); // CB
			refalpha_.push_back(unitconv_ * Vector3{ .0273,  .1680,  .2854}); // C
			refalpha_.push_back(unitconv_ * Vector3{ .0902,  .1993,  .3888}); // O
			refalpha_.push_back(unitconv_ * Vector3{ .0119,  .2532,  .1813}); // N    i+2
			refalpha_.push_back(unitconv_ * Vector3{ .0683,  .3916,  .1680}); // CA
			refalpha_.push_back(unitconv_ * Vector3{ .1580,  .3940,  .0395}); // CB
			refalpha_.push_back(unitconv_ * Vector3{-.0394,  .5011,  .1630}); // C
			refalpha_.push_back(unitconv_ * Vector3{-.1459,  .4814,  .0982}); // O
			refalpha_.push_back(unitconv_ * Vector3{-.2962,  .3559, -.1359}); // N    j - 2
			refalpha_.push_back(unitconv_ * Vector3{-.2439,  .2526, -.2287}); // CA
			refalpha_.push_back(unitconv_ * Vector3{-.1189,  .3006, -.3087}); // CB
			refalpha_.push_back(unitconv_ * Vector3{-.2081,  .1231, -.1520}); // C
			refalpha_.push_back(unitconv_ * Vector3{-.1524,  .1324, -.0409}); // O
			refalpha_.push_back(unitconv_ * Vector3{-.2326,  .0037, -.2095}); // N    j - 1
			refalpha_.push_back(unitconv_ * Vector3{-.1858, -.1269, -.1554}); // CA
			refalpha_.push_back(unitconv_ * Vector3{-.3053, -.2199, -.1291}); // CB
			refalpha_.push_back(unitconv_ * Vector3{-.0869, -.1949, -.2512}); // C
			refalpha_.push_back(unitconv_ * Vector3{-.1255, -.2070, -.3710}); // O
			refalpha_.push_back(unitconv_ * Vector3{ .0326, -.2363, -.2072}); // N    j
			refalpha_.push_back(unitconv_ * Vector3{ .1405, -.2992, -.2872}); // CA
			refalpha_.push_back(unitconv_ * Vector3{ .2699, -.2129, -.2917}); // CB
			refalpha_.push_back(unitconv_ * Vector3{ .1745, -.4399, -.2330}); // C
			refalpha_.push_back(unitconv_ * Vector3{ .1899, -.4545, -.1102}); // O
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

		static AntiBetaRMSDCV* Construct(const Json::Value& json, const std::string& path)
		{
			Json::ObjectRequirement validator;
			Json::Value schema;
			Json::Reader reader;

			reader.parse(JsonSchema::AntiBetaRMSDCV, schema);
			validator.Parse(schema, path);

			//Validate inputs
			validator.Validate(json, path);
			if(validator.HasErrors())
					throw BuildException(validator.GetErrors());

			std::vector<int> resids;
			for(auto& s : json["residue_ids"])
				resids.push_back(s.asInt());
			auto reference = json["reference"].asString();

			double unitconv = json.get("length_unit", 1).asDouble();

			int mode = json.get("mode", 0).asInt();

			return new AntiBetaRMSDCV(resids, reference, unitconv, mode);
		}
		//! Serialize this CV for restart purposes.
		/*!
		 * \param json JSON value
		 */
		virtual void Serialize(Json::Value& json) const override
		{
			json["type"] = "AntiBetaRMSD";
			json["reference"] = refpdb_;
			for(size_t i=0; i < resids_.size(); ++i)
				json["residue_ids"].append(resids_[i]);
			json["length_unit"] = unitconv_;
		}
	};
}

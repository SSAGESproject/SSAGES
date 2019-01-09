/**
 * This file is part of
 * SSAGES - Software Suite for Advanced General Ensemble Simulations
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

		//! Residue IDs for secondary structure calculation
		std::vector<int> resids_;

		//! Atom IDs for secondary structure calculation: backbone of resids_
		std::vector< std::vector<std::string> > atomids_;

		//! Atom IDs only, for secondary structure calculation: backbone atoms of resids_
		std::vector<int> ids_only_;

		//! Name of pdb reference for system
		std::string refpdb_;

		//! Length unit conversion: convert 1 nm to your internal MD units (ex. if using angstroms use 10)
		double unitconv_;

		//! Pairwise distance between backbone atoms for reference structure.
		std::vector<double> refdists_;

	public:
		//! Constructor.
		/*!
		 * \param resids IDs of residues for calculating secondary structure
		 * \param refpdb String of pdb filename with atom and residue indices.
		 * \param unitconv Conversion for internal MD length unit: 1 nm is equal to unitconv internal units
		 *
		 * Construct an AlphaRMSD CV -- calculates alpha helix character by
		 * summing pairwise RMSD to an ideal helix structure for all possible
		 * 6 residue segments.
		 *
		 */
		AlphaRMSDCV(std::vector<int> resids, std::string refpdb, double unitconv) :
		resids_(resids), refpdb_(refpdb), unitconv_(unitconv)
		{
			if(resids_.size() != 2 ){
				throw std::invalid_argument("AlphaRMSDCV: Input must designate range of residues with 2 residue numbers.");
			}

			resids_.clear();

			if(resids[0] >= resids[1]){
				throw std::invalid_argument("AlphaRMSDCV: Input must list lower residue index first: please reverse residue range.");
			} else if(resids[1] - resids[0] < 5) {
				throw std::invalid_argument("AlphaRMSDCV: Residue range must span at least 6 residues for alpha helix calculation.");
			}

			std::cout << "AlphaRMSDCV: Calculating alpha helix character from residue " << resids[0] << " to " << resids[1]  << "." << std::endl;

			for(int i = resids[0]; i <= resids[1]; i++){
				resids_.push_back(i);
			}
		}

		//! Initialize necessary variables.
		/*!
		 * \param snapshot Current simulation snapshot.
		 */
		void Initialize(const Snapshot& snapshot) override
		{
			atomids_ = ReadBackbone::GetPdbBackbone(refpdb_, resids_);
			ids_only_.resize(atomids_[0].size());

			// check to find all necessary atoms
			std::transform(atomids_[0].begin(), atomids_[0].end(), ids_only_.begin(), [](std::string val) {return std::stoi(val);});
			using std::to_string;
			std::vector<int> found;
			snapshot.GetLocalIndices(ids_only_, &found);
			size_t nfound = found.size();
			MPI_Allreduce(MPI_IN_PLACE, &nfound, 1, MPI_INT, MPI_SUM, snapshot.GetCommunicator());
			if(nfound != resids_.size() * 5)
				throw BuildException({
					"AlphaRMSDCV: Expected to find " +
					to_string(resids_.size() * 5) +
					" atoms, but only found " +
					to_string(nfound) + "."
				});

			// reference 'ideal' alpha helix structure, in nanometers
			// Reference structure taken from values used in Plumed 2.2.0
			std::vector<Vector3> refalpha;
			refalpha.push_back(unitconv_ * Vector3{ .0733,  .0519,  .5298 }); // N
			refalpha.push_back(unitconv_ * Vector3{ .1763,  .0810,  .4301 }); // CA
			refalpha.push_back(unitconv_ * Vector3{ .3166,  .0543,  .4881 }); // CB
			refalpha.push_back(unitconv_ * Vector3{ .1527, -.0045,  .3053 }); // C
			refalpha.push_back(unitconv_ * Vector3{ .1646,  .0436,  .1928 }); // O
			refalpha.push_back(unitconv_ * Vector3{ .1180, -.1312,  .3254 }); // N
			refalpha.push_back(unitconv_ * Vector3{ .0924, -.2203,  .2126 }); // CA
			refalpha.push_back(unitconv_ * Vector3{ .0650, -.3626,  .2626 }); // CB
			refalpha.push_back(unitconv_ * Vector3{-.0239, -.1711,  .1261 }); // C
			refalpha.push_back(unitconv_ * Vector3{-.0190, -.1815,  .0032 }); // O
			refalpha.push_back(unitconv_ * Vector3{-.1280, -.1172,  .1891 }); // N
			refalpha.push_back(unitconv_ * Vector3{-.2416, -.0661,  .1127 }); // CA
			refalpha.push_back(unitconv_ * Vector3{-.3548, -.0217,  .2056 }); // CB
			refalpha.push_back(unitconv_ * Vector3{-.1964,  .0529,  .0276 }); // C
			refalpha.push_back(unitconv_ * Vector3{-.2364,  .0659, -.0880 }); // O
			refalpha.push_back(unitconv_ * Vector3{-.1130,  .1391,  .0856 }); // N
			refalpha.push_back(unitconv_ * Vector3{-.0620,  .2565,  .0148 }); // CA
			refalpha.push_back(unitconv_ * Vector3{ .0228,  .3439,  .1077 }); // CB
			refalpha.push_back(unitconv_ * Vector3{ .0231,  .2129, -.1032 }); // C
			refalpha.push_back(unitconv_ * Vector3{ .0179,  .2733, -.2099 }); // O
			refalpha.push_back(unitconv_ * Vector3{ .1028,  .1084, -.0833 }); // N
			refalpha.push_back(unitconv_ * Vector3{ .1872,  .0593, -.1919 }); // CA
			refalpha.push_back(unitconv_ * Vector3{ .2850, -.0462, -.1397 }); // CB
			refalpha.push_back(unitconv_ * Vector3{ .1020,  .0020, -.3049 }); // C
			refalpha.push_back(unitconv_ * Vector3{ .1317,  .0227, -.4224 }); // O
			refalpha.push_back(unitconv_ * Vector3{-.0051, -.0684, -.2696 }); // N
			refalpha.push_back(unitconv_ * Vector3{-.0927, -.1261, -.3713 }); // CA
			refalpha.push_back(unitconv_ * Vector3{-.1933, -.2219, -.3074 }); // CB
			refalpha.push_back(unitconv_ * Vector3{-.1663, -.0171, -.4475 }); // C
			refalpha.push_back(unitconv_ * Vector3{-.1916, -.0296, -.5673 }); // O

			// calculate all necessary pairwise distances for reference structure
			std::fill(refdists_.begin(), refdists_.end(), 0.0);
			refdists_.resize(900, 0.0);
			for(size_t k = 0; k < 29; k++){
				for(size_t h = k + 1; h < 30; h++){
					refdists_[29 * k + h] = (refalpha[k] - refalpha[h]).norm();
				}
			}
		}

		//! Evaluate the CV.
		/*!
		 * \param snapshot Current simulation snapshot.
		 */
		void Evaluate(const Snapshot& snapshot) override
		{
			// need atom positions for all atoms
			const auto& pos = snapshot.GetPositions();
			auto& comm = snapshot.GetCommunicator();

			double rmsd, dist_norm, dxgrouprmsd;
			int localidx;
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
					localidx = snapshot.GetLocalIndex(ids_only_[5 * i + j]);
					if(localidx != -1) refxyz[j] = pos[localidx];
				}

				MPI_Allreduce(MPI_IN_PLACE, refxyz.data(), 90, MPI_DOUBLE, MPI_SUM, comm);

				// sum over 435 pairs in refxyz and refalpha_ to calculate CV
				for(size_t j = 0; j < 29; j++){
					for(size_t k = j + 1; k < 30; k++){
						dist_xyz = refxyz[j] - refxyz[k];
						dist_norm = dist_xyz.norm() - refdists_[29 * j + k];
						rmsd += dist_norm * dist_norm;
						deriv[j][k] = dist_xyz * dist_norm / dist_xyz.norm();
					}
				}

				rmsd = pow(rmsd / 435, 0.5) / 0.1;
				val_ += (1 - pow(rmsd, 8.0)) / (1 - pow(rmsd, 12.0));

				dxgrouprmsd = pow(rmsd, 11.0) + pow(rmsd, 7.0);
				dxgrouprmsd /= pow(rmsd, 8.0) + pow(rmsd, 4.0) + 1;
				dxgrouprmsd /= pow(rmsd, 8.0) + pow(rmsd, 4.0) + 1;
				dxgrouprmsd *= -40. / 435;

				for(size_t j = 0; j < 29; j++){
					for(size_t k = j + 1; k < 30; k++){
						localidx = snapshot.GetLocalIndex(ids_only_[5 * i + j]);
						if(localidx != -1) grad_[localidx] += dxgrouprmsd * deriv[j][k];
						localidx = snapshot.GetLocalIndex(ids_only_[5 * i + k]);
						if(localidx != -1) grad_[localidx] += dxgrouprmsd * deriv[j][k];
					}
				}
			}
		}

		//! \copydoc CollectiveVariable::BuildCV()
		static AlphaRMSDCV* Build(const Json::Value& json, const std::string& path)
		{
			Json::ObjectRequirement validator;
			Json::Value schema;
			Json::CharReaderBuilder rbuilder;
			Json::CharReader* reader = rbuilder.newCharReader();

			reader->parse(JsonSchema::AlphaRMSDCV.c_str(),
			              JsonSchema::AlphaRMSDCV.c_str() + JsonSchema::AlphaRMSDCV.size(),
			              &schema, nullptr);
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

			return new AlphaRMSDCV(resids, reference, unitconv);
		}
	};
}

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
		std::vector<int> _resids;

		//!< Atom IDs for secondary structure calculation: backbone of _resids
		std::vector<int> _atomids;

		//!< Name of pdb reference for system
		std::string _refpdb;		

		//!< Coordinates for reference structure
		std::vector<Vector3> _refalpha;

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
		AlphaRMSDCV(std::vector<int> resids, std::string refpdb) :
		_resids(resids), _refpdb(refpdb)
		{
			if(_resids.size() != 2 ){
				std::cout << "AlphaRMSDCV: Input must designate range of residues with 2 residue numbers." << std::endl;
				exit(0);
			}

			_resids.clear();

			if(resids[0] >= resids[1]){
				std::cout << "AlphaRMSDCV: Input must list lower residue index first: please reverse residue range." << std::endl;
				exit(0);
			} else if(resids[1] - resids[0] < 6) {
				std::cout << "AlphaRMSDCV: Residue range must span at least 6 residues for alpha helix calculation." << std::endl;
				exit(0);
			}
			
			std::cout << "AlphaRMSDCV: Calculating alpha helix character from residue " << resids[0] << " to " << resids[1]  << "." << std::endl;
			
			for(unsigned int i = resids[0]; i <= resids[1]; i++){
				_resids.push_back(i);
			}
			
		}

		// Initialize variables
		void Initialize(const Snapshot& snapshot) override
		{
			_atomids = ReadBackbone::GetPdbBackbone(_refpdb, _resids);

			// Reference 'ideal' alpha helix structure -- in angstroms
			//_refalpha.push_back(Vector3{ 0.733,  0.519,  5.298 }); // N    
			//_refalpha.push_back(Vector3{ 1.763,  0.810,  4.301 }); // CA
			//_refalpha.push_back(Vector3{ 3.166,  0.543,  4.881 }); // CB
			//_refalpha.push_back(Vector3{ 1.527, -0.045,  3.053 }); // C
			//_refalpha.push_back(Vector3{ 1.646,  0.436,  1.928 }); // O
			//_refalpha.push_back(Vector3{ 1.180, -1.312,  3.254 }); // N    
			//_refalpha.push_back(Vector3{ 0.924, -2.203,  2.126 }); // CA
			//_refalpha.push_back(Vector3{ 0.650, -3.626,  2.626 }); // CB
			//_refalpha.push_back(Vector3{-0.239, -1.711,  1.261 }); // C
			//_refalpha.push_back(Vector3{-0.190, -1.815,  0.032 }); // O
			//_refalpha.push_back(Vector3{-1.280, -1.172,  1.891 }); // N   
			//_refalpha.push_back(Vector3{-2.416, -0.661,  1.127 }); // CA
			//_refalpha.push_back(Vector3{-3.548, -0.217,  2.056 }); // CB
			//_refalpha.push_back(Vector3{-1.964,  0.529,  0.276 }); // C
			//_refalpha.push_back(Vector3{-2.364,  0.659, -0.880 }); // O
			//_refalpha.push_back(Vector3{-1.130,  1.391,  0.856 }); // N   
			//_refalpha.push_back(Vector3{-0.620,  2.565,  0.148 }); // CA
			//_refalpha.push_back(Vector3{ 0.228,  3.439,  1.077 }); // CB
			//_refalpha.push_back(Vector3{ 0.231,  2.129, -1.032 }); // C
			//_refalpha.push_back(Vector3{ 0.179,  2.733, -2.099 }); // O
			//_refalpha.push_back(Vector3{ 1.028,  1.084, -0.833 }); // N   
			//_refalpha.push_back(Vector3{ 1.872,  0.593, -1.919 }); // CA
			//_refalpha.push_back(Vector3{ 2.850, -0.462, -1.397 }); // CB
			//_refalpha.push_back(Vector3{ 1.020,  0.020, -3.049 }); // C
			//_refalpha.push_back(Vector3{ 1.317,  0.227, -4.224 }); // O
			//_refalpha.push_back(Vector3{-0.051, -0.684, -2.696 }); // N   
			//_refalpha.push_back(Vector3{-0.927, -1.261, -3.713 }); // CA
			//_refalpha.push_back(Vector3{-1.933, -2.219, -3.074 }); // CB
			//_refalpha.push_back(Vector3{-1.663, -0.171, -4.475 }); // C
			//_refalpha.push_back(Vector3{-1.916, -0.296, -5.673 }); // O
			
			// in nanometers
			_refalpha.push_back(Vector3{ .0733,  .0519,  .5298 }); // N    
			_refalpha.push_back(Vector3{ .1763,  .0810,  .4301 }); // CA
			_refalpha.push_back(Vector3{ .3166,  .0543,  .4881 }); // CB
			_refalpha.push_back(Vector3{ .1527, -.0045,  .3053 }); // C
			_refalpha.push_back(Vector3{ .1646,  .0436,  .1928 }); // O
			_refalpha.push_back(Vector3{ .1180, -.1312,  .3254 }); // N    
			_refalpha.push_back(Vector3{ .0924, -.2203,  .2126 }); // CA
			_refalpha.push_back(Vector3{ .0650, -.3626,  .2626 }); // CB
			_refalpha.push_back(Vector3{-.0239, -.1711,  .1261 }); // C
			_refalpha.push_back(Vector3{-.0190, -.1815,  .0032 }); // O
			_refalpha.push_back(Vector3{-.1280, -.1172,  .1891 }); // N   
			_refalpha.push_back(Vector3{-.2416, -.0661,  .1127 }); // CA
			_refalpha.push_back(Vector3{-.3548, -.0217,  .2056 }); // CB
			_refalpha.push_back(Vector3{-.1964,  .0529,  .0276 }); // C
			_refalpha.push_back(Vector3{-.2364,  .0659, -.0880 }); // O
			_refalpha.push_back(Vector3{-.1130,  .1391,  .0856 }); // N   
			_refalpha.push_back(Vector3{-.0620,  .2565,  .0148 }); // CA
			_refalpha.push_back(Vector3{ .0228,  .3439,  .1077 }); // CB
			_refalpha.push_back(Vector3{ .0231,  .2129, -.1032 }); // C
			_refalpha.push_back(Vector3{ .0179,  .2733, -.2099 }); // O
			_refalpha.push_back(Vector3{ .1028,  .1084, -.0833 }); // N   
			_refalpha.push_back(Vector3{ .1872,  .0593, -.1919 }); // CA
			_refalpha.push_back(Vector3{ .2850, -.0462, -.1397 }); // CB
			_refalpha.push_back(Vector3{ .1020,  .0020, -.3049 }); // C
			_refalpha.push_back(Vector3{ .1317,  .0227, -.4224 }); // O
			_refalpha.push_back(Vector3{-.0051, -.0684, -.2696 }); // N   
			_refalpha.push_back(Vector3{-.0927, -.1261, -.3713 }); // CA
			_refalpha.push_back(Vector3{-.1933, -.2219, -.3074 }); // CB
			_refalpha.push_back(Vector3{-.1663, -.0171, -.4475 }); // C
			_refalpha.push_back(Vector3{-.1916, -.0296, -.5673 }); // O
		}

		// Evaluate the CV
		void Evaluate(const Snapshot& snapshot) override
		{
			// need atom positions for all atoms in _atomids
			const auto& pos = snapshot.GetPositions();
			//const auto& image_flags = snapshot.GetImageFlags();
			std::vector<int> groupidx;
			snapshot.GetLocalIndices(_atomids, &groupidx);	// get correct local atom indices
			
			unsigned int resgroups = _resids.size() - 5;

			double rmsd, dist_norm, dxgrouprmsd;
			Vector3 dist_xyz, dist_ref;
			std::vector<Vector3> _refxyz;
			std::vector< std::vector< Vector3 > > deriv(30, std::vector<Vector3>(30, Vector3{0,0,0})); 
			
			std::fill(_grad.begin(), _grad.end(), Vector3{0,0,0});
			_grad.resize(snapshot.GetNumAtoms(), Vector3{0,0,0});
			_val = 0.0;

			// for each set of 6 residues
			for(size_t i = 0; i < resgroups; i++){

				// clear temp rmsd calculation
				rmsd = 0.0;
				
				// load _refxyz with the correct 30 reference atoms
				std::fill(_refxyz.begin(), _refxyz.end(), Vector3{0,0,0});
				_refxyz.resize(30, Vector3{0,0,0});
				for(size_t j = 0; j < 30; j++){
					_refxyz[j] = pos[groupidx[5 * i + j]];
				}

				// sum over 435 pairs in _refxyz and _refalpha to calculate CV
				for(size_t j = 0; j < 29; j++){
					for(size_t k = j + 1; k < 30; k++){
						dist_xyz = _refxyz[j] - _refxyz[k];
						dist_ref = _refalpha[j] - _refalpha[k]; // could be precalculated
						dist_norm = dist_xyz.norm() - dist_ref.norm();
						rmsd += dist_norm * dist_norm;
						
						deriv[j][k] = dist_xyz * dist_norm / dist_xyz.norm();
					}
				}

				rmsd /= 43.5;
				_val += (1 - pow(rmsd, 8.0)) / (1 - pow(rmsd, 12.0));

				dxgrouprmsd = 8E8 * (5E3 * pow(rmsd, 11.0) + pow(rmsd, 7.0));
				dxgrouprmsd /= 1E8 * pow(rmsd, 8.0) + 1E4 * pow(rmsd, 4.0) + 1;
				dxgrouprmsd /= 1E8 * pow(rmsd, 8.0) + 1E4 * pow(rmsd, 4.0) + 1;
				dxgrouprmsd *= -2.0 / 435;  // derivative has prefactor of -2 / (30 choose 2)

				for(size_t j = 0; j < 29; j++){
					for(size_t k = j + 1; k < 30; k++){
						_grad[groupidx[5 * i + j]] += dxgrouprmsd * deriv[j][k];
						_grad[groupidx[5 * i + k]] -= dxgrouprmsd * deriv[j][k];
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
			json["reference"] = _refpdb;
			for(size_t i=0; i < _resids.size(); ++i)
				json["residue_ids"].append(_resids[i]);
		}
	};
}

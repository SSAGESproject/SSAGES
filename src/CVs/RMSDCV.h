/**
 * This file is part of
 * SSAGES - Suite for Advanced Generalized Ensemble Simulations
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

#include "CollectiveVariable.h"
#include "../Utility/ReadFile.h"
#include <array>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

namespace SSAGES
{
	//! Collective variable to calculate root mean square displacement.
	/*!
	 *RMSD calculation performed as reported in 
	 *Coutsias, E. A., Seok, C., and Dill, K. A., "Using Quaternions to Calculate RMSD", J. Comput. Chem. 25: 1849-1857, 2004
	 */

	class RMSDCV : public CollectiveVariable
	{
	private:

		std::vector<int> _atomids; //!< IDs of the atoms used for Rg calculation

		std::vector<int> _pertatoms; //!< Array to store indicies of atoms of interest

		//! Name of model structure.
		std::string _molecule; 

		//! Store reference structure coordinates.
		std::vector<Vector3> _refcoord;

		//! Center of mass of reference.
		Vector3 _COMref;

	public:
		//! Constructor.
		/*!
		 * \param atomids IDs of the atoms defining Rg.
		 * \param molxyz String determining the molecule.
		 * \param use_range If \c True Use range of atoms defined by the two atoms in atomids.
		 *
		 * Construct a RMSD CV.
		 *
		 * \todo Bounds needs to be an input and periodic boundary conditions
		 */
		RMSDCV(std::vector<int> atomids, std::string molxyz, bool use_range = false) :
		_atomids(atomids), _molecule(molxyz)
		{
			if(use_range)
			{
				if(atomids.size() != 2)
				{	std::cout<<"RMSDCV: If using range, must define only two atoms!"<<std::endl;
					exit(0);
				}

				_atomids.clear();

				if(atomids[0] >= atomids[1])
				{	std::cout<<"RMSDCV: Please reverse atom range or check that atom range is not equal!"<<std::endl;
					exit(0);
				}					
				for(int i = atomids[0]; i <= atomids[1];i++)
					_atomids.push_back(i);
			}
		}

		// Initialize variables
		void Initialize(const Snapshot& snapshot) override
		{

			const auto& mass = snapshot.GetMasses();
			const auto& ids = snapshot.GetAtomIDs();

			// Initialize gradient
			auto n = snapshot.GetPositions().size();
			_grad.resize(n);
			_refcoord.resize(_atomids.size());

			std::vector<std::array<double,4>> xyzinfo = ReadFile::ReadXYZ(_molecule);

			if(_refcoord.size() != xyzinfo.size())
				throw std::runtime_error("Reference structure and input structure atom size do not match!");

			for(size_t i = 0; i < xyzinfo.size(); i++)
			{
				_refcoord[i][0] = xyzinfo[i][1];
				_refcoord[i][1] = xyzinfo[i][2];
				_refcoord[i][2] = xyzinfo[i][3];
			}

			Vector3 mass_pos_prod_ref;
			mass_pos_prod_ref.setZero();
			double total_mass = 0;

			// Loop through atom positions
			for( size_t i = 0; i < mass.size(); ++i)
			{
				// Loop through pertinent atoms
				for(size_t j = 0; j < _atomids.size(); j++)
				{
					if(ids[i] == _atomids[j])
					{
						mass_pos_prod_ref += mass[i]*_refcoord[j];
						total_mass += mass[i];
						break;
					}
				}
			}

			_COMref = mass_pos_prod_ref/total_mass;
		}

		// Evaluate the CV
		void Evaluate(const Snapshot& snapshot) override
		{
			const auto& pos = snapshot.GetPositions();
			const auto& ids = snapshot.GetAtomIDs();
			const auto& mass = snapshot.GetMasses();
			const auto& image_flags = snapshot.GetImageFlags();
			Vector3 mass_pos_prod = {0,0,0};
			double total_mass=0;
			int i;
			// Loop through atom positions
			for( size_t i = 0; i < pos.size(); ++i)
			{
				_grad[i].setZero();
				// Loop through pertinent atoms
				for(size_t j = 0; j < _atomids.size(); j++)
				{
					if(ids[i] == _atomids[j])
					{
						_pertatoms[j] = i;
						auto u_coord = snapshot.UnwrapVector(pos[i], image_flags[i]);

						mass_pos_prod += mass[i]*u_coord;
						total_mass += mass[i];
						break;
					}
				}
			}


			// Calculate center of mass
			//Vector3 mass_pos_prod;
			//mass_pos_prod.setZero();
			total_mass = 0;
			Vector3 _COM;

			// Need to unwrap coordinates for appropriate COM calculation. 

			//for( size_t j = 0; j < _pertatoms.size(); ++j)
			//{
			//	i = _pertatoms[j];
			//	auto u_coord = UnwrapCoordinates(snapshot.GetLatticeConstants(), image_flags[i], pos[i]);

			//	mass_pos_prod[0] += mass[i]*u_coord[0];
			//	mass_pos_prod[1] += mass[i]*u_coord[1];
			//	mass_pos_prod[2] += mass[i]*u_coord[2];
			//	total_mass += mass[i];
			//}

			_COM = mass_pos_prod/total_mass;

			// Build correlation matrix
			Matrix3 R;

			Vector3 diff;
			Vector3 diff_ref;
			double part_RMSD = 0;

			for( size_t j = 0; j < _pertatoms.size(); ++j)
			{
				i = _pertatoms[j];
				auto u_coord = snapshot.UnwrapVector(pos[i], image_flags[i]);

				diff = u_coord -_COM;
				diff_ref = _refcoord[j] - _COMref;

				R(0,0) += diff[0]*diff_ref[0];
				R(0,1) += diff[0]*diff_ref[1];
				R(0,2) += diff[0]*diff_ref[2];
				R(1,0) += diff[1]*diff_ref[0];
				R(1,1) += diff[1]*diff_ref[1];
				R(1,2) += diff[1]*diff_ref[2];
				R(2,0) += diff[2]*diff_ref[0];
				R(2,1) += diff[2]*diff_ref[1];
				R(2,2) += diff[2]*diff_ref[2];


				auto normdiff2 = diff.norm()*diff.norm();
				auto normref2 = diff_ref.norm()*diff_ref.norm();

				part_RMSD += (normdiff2 + normref2);
			}

			Eigen::Matrix4d F;
			F(0,0)= R(0,0) + R(1,1) + R(2,2);
			F(1,0)= R(1,2) - R(2,1);
			F(2,0)= R(2,0) - R(0,2);
			F(3,0)= R(0,1) - R(1,0);
			
			F(0,1)= R(1,2) - R(2,1);
			F(1,1)= R(0,0) - R(1,1) - R(2,2);
			F(2,1)= R(0,1) + R(1,0);
			F(3,1)= R(0,2) + R(2,0);

			F(0,2)= R(2,0) - R(0,2);
			F(1,2)= R(0,1) - R(1,0);
			F(2,2)= -R(0,0) + R(1,1) - R(2,2);
			F(3,2)= R(1,2) + R(2,1);

			F(0,3)= R(0,1) - R(1,0);
			F(1,3)= R(0,2) - R(2,0);
			F(2,3)= R(1,2) - R(2,1);
			F(3,3)= -R(0,0) - R(1,1) + R(2,2);
			
			//Find eigenvalues
			Eigen::EigenSolver<Eigen::Matrix4d> es(F);
			//EigenSolver<F> es;
			//es.compute(F);
			//auto max_lambda = es.eigenvalues().maxCoeff();
			auto lambda = es.eigenvalues().real();
			//double max_lambda;
			auto max_lambda = lambda.maxCoeff();
			int column;
			

			for (i =0; i < 4; ++i)
				if(es.eigenvalues()[i] == max_lambda)
					column=i;

			// Calculate RMSD


			auto RMSD = sqrt((part_RMSD - (2*max_lambda))/(_atomids.size()));

			_val = RMSD;

			// Calculate gradient

			auto eigenvector = es.eigenvectors().col(column);

			auto q0 = eigenvector[0].real();
			auto q1 = eigenvector[1].real();
			auto q2 = eigenvector[2].real();
			auto q3 = eigenvector[3].real();

			Matrix3 RotMatrix(3,3);
			RotMatrix(0,0) = q0*q0+q1*q1-q2*q2-q3*q3;
			RotMatrix(1,0) = 2*(q1*q2+q0*q3);
			RotMatrix(2,0) = 2*(q1*q3-q0*q2);

			RotMatrix(0,1) = 2*(q1*q2-q0*q3);
			RotMatrix(1,1) = q0*q0-q1*q1+q2*q2-q3*q3;
			RotMatrix(2,1) = 1*(q2*q3+q0*q1);

			RotMatrix(0,2) = 2*(q1*q3+q0*q2);
			RotMatrix(1,2) = 2*(q2*q3-q0*q1);
			RotMatrix(2,2) = q0*q0-q1*q1-q2*q2+q3*q3;

			//double trans[3];
			for(size_t j=0; j<_pertatoms.size();++j)
			{
				auto trans = RotMatrix.transpose()*_refcoord[j];

				i = _pertatoms[j];
				auto u_coord = pos[i]; //UnwrapCoordinates(snapshot.GetLatticeConstants(), image_flags[i], pos[i]);

				_grad[i][0] = (u_coord[0] - trans[0])/((_atomids.size())*_val);
				_grad[i][1] = (u_coord[1] - trans[1])/((_atomids.size())*_val);
				_grad[i][2] = (u_coord[2] - trans[2])/((_atomids.size())*_val);
			}
		}

		//! Serialize this CV for restart purposes.
		/*!
		 * \param json JSON value
		 */
		virtual void Serialize(Json::Value& json) const override
		{
			json["type"] = "RMSD";
			json["reference"] = _molecule;
			for(size_t i=0; i < _atomids.size(); ++i)
				json["atom ids"].append(_atomids[i]);
			for(size_t i = 0; i < _bounds.size(); ++i)
				json["bounds"].append(_bounds[i]);
		}
	};
}

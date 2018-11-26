/**
 * This file is part of
 * SSAGES - Software Suite for Advanced General Ensemble Simulations
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
#include "Validator/ObjectRequirement.h"
#include "Drivers/DriverException.h"
#include "Snapshot.h"
#include "schema.h"
#include <Eigen/Eigenvalues>

namespace SSAGES
{
	//! Define components of gyration tensor
	/*!
	 * Components of the gyration tensor. CV uses these definitions
	 * to determine what to calculate.
	 */
	enum GyrationTensor
	{
		Rg = 0, // Radius of gyration (squared)
		principal1 = 1, // First (largest) principal moment 
		principal2 = 2, // Second (middle) principal moment 
		principal3 = 3, // Third (smallest) principal moment
		asphericity = 4, // Asphericity 
		acylindricity = 5, // Acylindricity
		shapeaniso = 6 // Relative shape anisotropy
	};

	//! Collective variable on components of the gyration tensor.
	/*!
	 * Collective variable on components of gyration tensor. Depending on the 
	 * user selection, this will specify a principal moment, radius of gyration, 
	 * or another shape descriptor.
	 *
	 * \ingroup CVs
	 */
	class GyrationTensorCV : public CollectiveVariable
	{
	private:
		Label atomids_; //!< IDs of the atoms used for calculation
		GyrationTensor component_; //!< Component of gyration tensor to compute.

		//! Each dimension determines if it is applied by the CV.
		Bool3 dim_;

	public:
		//! Constructor.
		/*!
		 * \param atomids IDs of the atoms defining gyration tensor.
		 * \param component Specification of component to compute.
		 *
		 * Construct a GyrationTensorCV.
		 *
		 */
		GyrationTensorCV(const Label& atomids, GyrationTensor component) :
		atomids_(atomids), component_(component), dim_{true, true, true}
		{
		}

		//! Constructor.
		/*!
		 * \param atomids IDs of the atoms defining gyration tensor.
		 * \param component Specification of component to compute.
		 * \param dimx If \c True, include x dimension.
		 * \param dimy If \c True, include y dimension.
		 * \param dimz If \c True, include z dimension.
		 *
		 * Construct a GyrationTensorCV.
		 *
		 */
		GyrationTensorCV(const Label& atomids, GyrationTensor component, bool dimx, bool dimy, bool dimz) :
		atomids_(atomids), component_(component), dim_{dimx, dimy, dimz}
		{
		}

		//! Initialize necessary variables.
		/*!
		 * \param snapshot Current simulation snapshot.
		 */
		void Initialize(const Snapshot& snapshot) override
		{
			using std::to_string;

			auto n = atomids_.size();

			// Make sure atom ID's are on at least one processor. 
			std::vector<int> found(n, 0);
			for(size_t i = 0; i < n; ++i)
			{
				if(snapshot.GetLocalIndex(atomids_[i]) != -1)
					found[i] = 1;
			}

			MPI_Allreduce(MPI_IN_PLACE, found.data(), n, MPI_INT, MPI_SUM, snapshot.GetCommunicator());
			unsigned ntot = std::accumulate(found.begin(), found.end(), 0, std::plus<int>());
			if(ntot != n)
				throw BuildException({
					"GyrationTensorCV: Expected to find " + 
					to_string(n) + 
					" atoms, but only found " + 
					to_string(ntot) + "."
				});		
		}		

		//! Evaluate the CV.
		/*!
		 * \param snapshot Current simulation snapshot.
		 */
		void Evaluate(const Snapshot& snapshot) override
		{
			using namespace Eigen; 

			// Get local atom indices and compute COM. 
			std::vector<int> idx;
			snapshot.GetLocalIndices(atomids_, &idx);

			// Get data from snapshot. 
			auto n = snapshot.GetNumAtoms();
			const auto& masses = snapshot.GetMasses();
			const auto& pos = snapshot.GetPositions();

			// Initialize gradient.
			std::fill(grad_.begin(), grad_.end(), Vector3{0,0,0});
			grad_.resize(n, Vector3{0,0,0});
			
			// Compute total and center of mass.
			auto masstot = snapshot.TotalMass(idx);
			Vector3 com = snapshot.CenterOfMass(idx, masstot);

			// Gyration tensor and temporary vector to store positions in inertial frame. 
			Matrix3 S = Matrix3::Zero();
			std::vector<Vector3> ris; 
			ris.reserve(idx.size());
			for(auto& i : idx)
			{
				ris.emplace_back(snapshot.ApplyMinimumImage(pos[i] - com).cwiseProduct(dim_.cast<double>()));
				S.noalias() += masses[i]*ris.back()*ris.back().transpose();
			}

			// Reduce gyration tensor across processors and normalize.
			MPI_Allreduce(MPI_IN_PLACE, S.data(), S.size(), MPI_DOUBLE, MPI_SUM, snapshot.GetCommunicator());
			S /= masstot;

			// Perform EVD. The columns are the eigenvectors. 
			// SelfAdjoint solver sorts in ascending order. 
			SelfAdjointEigenSolver<Matrix3> solver;
			solver.computeDirect(S);
			const auto& eigvals = solver.eigenvalues();
			const auto& eigvecs = solver.eigenvectors();
			
			// Assign variables for clarity. l1 is largest.
			auto l1 = eigvals[2], 
			     l2 = eigvals[1], 
			     l3 = eigvals[0];
			const auto& n1 = eigvecs.col(2), 
			            n2 = eigvecs.col(1), 
			            n3 = eigvecs.col(0);

			val_ = 0;
			auto sum = l1 + l2 + l3;
			auto sqsum = l1*l1 + l2*l2 + l3*l3;
			switch(component_)
			{
				case Rg:
					val_ = sum;
					break;
				case principal1:
					val_ = l1;
					break;
				case principal2:
					val_ = l2;
					break;
				case principal3:
					val_ = l3;
					break;
				case asphericity:
					val_ = l1 - 0.5*(l2 + l3);
					break;
				case acylindricity:
					val_ = l2 - l3;
					break;
				case shapeaniso:
					val_ = 1.5*sqsum/(sum*sum) - 0.5;
					break;
			}

			// Compute gradient.
			size_t j = 0;
			for(auto& i : idx)
			{
				// Compute derivative of each eigenvalue and use combos in components. 
				auto dl1 = 2.*masses[i]/masstot*(1.-masses[i]/masstot)*ris[j].dot(n1)*n1;
				auto dl2 = 2.*masses[i]/masstot*(1.-masses[i]/masstot)*ris[j].dot(n2)*n2;
				auto dl3 = 2.*masses[i]/masstot*(1.-masses[i]/masstot)*ris[j].dot(n3)*n3;

				switch(component_)
				{
					case Rg:
						grad_[i] = dl1 + dl2 + dl3;
						break;
					case principal1:
						grad_[i] = dl1;
						break;
					case principal2:
						grad_[i] = dl2;
						break;
					case principal3:
						grad_[i] = dl3;
						break;
					case asphericity:
						grad_[i] = dl1 - 0.5*(dl2 + dl3);
						break;
					case acylindricity:
						grad_[i] = dl2 - dl3;
						break;
					case shapeaniso:
						grad_[i] = 3.*(l1*dl1+l2*dl2+l3*dl3)/(sum*sum) - 
						           3.*sqsum*(dl1+dl2+dl3)/(sum*sum*sum);
						break;
				}

				++j;
			}
		}

		//! \copydoc CollectiveVariable::BuildCV()
		static GyrationTensorCV* Build(const Json::Value& json, const std::string& path)
		{
			Json::ObjectRequirement validator;
			Json::Value schema;
			Json::CharReaderBuilder rbuilder;
			Json::CharReader* reader = rbuilder.newCharReader();

			reader->parse(JsonSchema::GyrationTensorCV.c_str(),
			              JsonSchema::GyrationTensorCV.c_str() + JsonSchema::GyrationTensorCV.size(),
			              &schema, NULL);
			validator.Parse(schema, path);

			// Validate inputs. 
			validator.Validate(json, path); 
			if(validator.HasErrors())
				throw BuildException(validator.GetErrors());

			std::vector<int> atomids; 
			for(auto& s : json["atom_ids"])
				atomids.push_back(s.asInt());

			GyrationTensor component = Rg;
			auto comp = json["component"].asString();
			
			if(comp == "Rg")
				component = Rg; 
			else if(comp == "principal1")
				component = principal1;
			else if(comp == "principal2")
				component = principal2;
			else if(comp == "principal3")
				component = principal3;
			else if(comp == "asphericity")
				component = asphericity;
			else if(comp == "acylindricity")
				component = acylindricity;
			else if(comp == "shapeaniso")
				component = shapeaniso;

			GyrationTensorCV* c;
			if(json.isMember("dimension"))
			{
				auto dimx = json["dimension"][0].asBool();
				auto dimy = json["dimension"][1].asBool();
				auto dimz = json["dimension"][2].asBool();
				c = new GyrationTensorCV(atomids, component, dimx, dimy, dimz);
			}
			else
				c = new GyrationTensorCV(atomids, component);

			return c;
		}
	};
}

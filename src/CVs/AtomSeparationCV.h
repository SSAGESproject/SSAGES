/**
 * This file is part of
 * SSAGES - Suite for Advanced Generalized Ensemble Simulations
 *
 * Copyright 2016 Emre Sevgen <sesevgen@uchicago.edu>
 *                Hythem Sidky <hsidky@nd.edu>
 *                Ben Sikora <bsikora906@gmail.com>
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

#include <array>
#include <math.h>

namespace SSAGES
{
	//! Collective variable on the distance between two atoms.
	/*!
	 *  Collective variable on two atom positions. This CV will return the
	 *  distance between two specific atoms of the simulation.
	 *
	 *  \ingroup CVs
	 */
	class AtomSeparationCV : public CollectiveVariable
	{
	private:
		//! ID of atom of interest 1.
		int _atomid1;
		
		//! ID of atom of interest 1.
		int _atomid2;

		//! Current value of the CV.
		double _val;

		//! Gradient of the CV, dr/dxi.
		std::vector<Vector3> _grad;

		//! Bounds on CV.
		std::array<double, 2> _bounds;
		

	public:
		//! Constructor
		/*!
		 * \param atomid1 ID of the first atom of interest.
		 * \param atomid2 ID of the second atom of interest.
		 *
		 * Construct an atom separation CV.
		 *
		 * \todo bounds needs to be an input.
		 * \todo Consider non-periodic boundary conditions.
		 */
		AtomSeparationCV(int atomid1, int atomid2) : 
		_atomid1(atomid1),_atomid2(atomid2), _val(0), _grad(0), _bounds{{0,0}}
		{
		}

		//! Initialize necessary variables.
		/*!
		 * \param snapshot Current simulation snapshot.
		 */
		void Initialize(const Snapshot& snapshot) override
		{
			// Initialize gradient. 
			auto n = snapshot.GetPositions().size();		
			_grad.resize(n);
		}
		

		// Evaluate the CV.
		/*!
		 * \param snapshot Current simulation snapshot.
		 */
		void Evaluate(const Snapshot& snapshot) override
		{
			// Gradient and value. 
			const auto& pos = snapshot.GetPositions(); 
			const auto& ids = snapshot.GetAtomIDs();

			// Some temp variables.
			Vector3 pos1{0, 0, 0};
			Vector3 pos2{0, 0, 0};
			
			size_t index1 = 0, index2 = 0;		
			
			// Loop through atom positions
			for(size_t i = 0; i < pos.size(); ++i)
			{
				_grad[i].setZero();

				// If we are at the atom ID of interest.
				if(ids[i] == _atomid1)
				{ 
					pos1 = pos[i];
					index1 = i;
				}
				else if(ids[i] == _atomid2)
				{
					pos2 = pos[i];
					index2 = i;
				}
			}

			// Account for PBCs
			Vector3 rij = snapshot.ApplyMinimumImage(pos1 - pos2);	

			// Compute norm and gradient.
			auto r = rij.norm();
			_grad[index1] = rij/r;
			_grad[index2] = -rij/r;
			_val = r;
		}

		//! Return the value of the CV.
		/*!
		 * \return Current value of the CV.
		 */
		double GetValue() const override 
		{	 
			return _val; 
		}

		//! Return input value taking periodic boundary conditions into account.
		/*!
		 * \param Location Input value.
		 * \return Unmodified input value.
		 *
		 * As the distance has no periodic boundary conditions, the input value
		 * is not modified.
		 */
		double GetPeriodicValue(double Location) const override
		{
			return Location;
		}

		//! Return the gradient of the CV.
		/*!
		 * \return Gradient of the CV.
		 */
		const std::vector<Vector3>& GetGradient() const override
		{
			return _grad;
		}

		//! Return the boundaries of the CV.
		/*!
		 * \return Values of the lower and upper boundaries of the CV.
		 */
		const std::array<double, 2>& GetBoundaries() const override
		{
			return _bounds;
		}

		//! Get difference taking periodic boundary conditions into account.
		/*!
		 * \param Location Input value.
		 * \return Simple difference.
		 *
		 * Return the simple difference between the current value of the CV and
		 * the input value.
		 */
		double GetDifference(const double Location) const override
		{
			return _val - Location;
		}
		
		//! Serialize this CV for restart purposes.
		/*!
		 * \param json JSON value
		 */
		virtual void Serialize(Json::Value& json) const override
		{
			json["type"] = "AtomSeparation";
			json["atom id 1"] = _atomid1;
			json["atom id 2"] = _atomid2;
			for(size_t i = 0; i < _bounds.size(); ++i)
				json["bounds"].append(_bounds[i]);
		}
	};
}

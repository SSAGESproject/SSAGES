#pragma once 

#include "CollectiveVariable.h"

#include <array>
#include <cmath>

namespace SSAGES
{
	//! Collective variable on the torsion angles.
	/*!
	 * Collective variable on an proper dihedral. This will return the angle
	 * between two planes as defined in \cite VANSCHAIK1993751. Singularities
	 * are avoided as described in \cite BLONDEL19961132.
	 *
	 * \ingroup CVs
	 */

	class TorsionalCV : public CollectiveVariable
	{
	private:
		int _atomid1; //!< ID of the first atom defining the torsion angle.
		int _atomid2; //!< ID of the second atom defining the torsion angle.
		int _atomid3; //!< ID of the third atom defining the torsion angle.
		int _atomid4; //!< ID of the forth atom defining the torsion angle.

		//! Current value of the CV.
		double _val;

		//! If \c True, use periodic boundary conditions.
		bool _periodic;

		//! Gradients of the Dihedral CV, dtheta/dri, dtheta/drj, dtheta/drk, dtheta/drl.
		std::vector<Vector3> _grad;

		//! Bounds on CV.
		std::array<double, 2> _bounds;

	public:
		//! Constructor.
		/*!
		 * \param atomid1 ID of the first atom defining the dihedral angle.
		 * \param atomid2 ID of the second atom defining the dihedral angle.
		 * \param atomid3 ID of the third atom defining the dihedral angle.
		 * \param atomid4 ID of the forth atom defining the dihedral angle.
		 * \param periodic If \c True consider periodic boundary conditions.
		 *
		 * Construct an dihedral CV.
		 *
		 * \todo Bounds needs to be an input and periodic boundary conditions
		 */
		TorsionalCV(int atomid1, int atomid2, int atomid3, int atomid4, bool periodic) : 
		_atomid1(atomid1), _atomid2(atomid2), _atomid3(atomid3), _atomid4(atomid4),
		_val(0), _periodic(periodic), _grad(0), _bounds{{0,0}}
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

		//! Evaluate the CV.
		/*!
		 * \param snapshot Current simulation snapshot.
		 */
		void Evaluate(const Snapshot& snapshot) override
		{
			// Gradient and value. 
			const auto& pos = snapshot.GetPositions(); 
			const auto& ids = snapshot.GetAtomIDs();

			int iindex, jindex, kindex, lindex;

			iindex = jindex = kindex = lindex = -1;
			Vector3 xi, xj, xk, xl;
			
			// Loop through atom positions
			for(size_t i = 0; i < pos.size(); ++i)
			{
				_grad[i].setZero();
				// If we are at the atom ID of interest, grab coordinates
				if(ids[i] == _atomid1)
				{
					//coordinates for atom i
					xi = pos[i];
					iindex = i;
				}
				if(ids[i] == _atomid2)
				{
					//coordinates for atom j
					xj = pos[i];
					jindex = i;
				}
				if(ids[i] == _atomid3)
				{
					//coordinates for atom k
					xk = pos[i];
					kindex = i;
				}
				if(ids[i] == _atomid4)
				{
					//coordinates for atom l
					xl = pos[i];
					lindex = i;
				}
			}

			if(iindex < 0 || jindex < 0 || kindex < 0 || lindex <0)
			{
				std::cout<<"Out of bounds index, could not locate an ID"<<std::endl;
				exit(0);
			}

			//Calculate pertinent vectors
			auto F = xi - xj;
			auto G = xj - xk;
			auto H = xl - xk;
			auto A = F.cross(G);
			auto B = H.cross(G);
			
			auto y = B.cross(A).dot(G.normalized());
			auto x = A.dot(B);

			_val = atan2(y, x);


			auto Zed = F.dot(G.normalized())/A.dot(A); 
			auto Ned = H.dot(G.normalized())/B.dot(B);

			_grad[iindex] = -G.norm()*A/A.dot(A);
			_grad[lindex] = G.norm()*B/B.dot(B);
			_grad[jindex] = Zed*A - Ned*B - _grad[iindex];
			_grad[kindex] = Ned*B  - Zed*A - _grad[lindex];
		}

		//! Return the value of the CV.
		/*!
		 * \return Current value of the CV.
		 */
		double GetValue() const override 
		{ 
			return _val; 
		}

		//! Return value taking periodic boundary conditions into account
		/*!
		 * \param Location Input value.
		 * \return Wrapped or unwrapped input value depending on whether
		 *         periodic boundaries are used.
		 *
		 * If periodic boundaries are used, this function wraps the input
		 * value into the range (-pi, pi). Otherwise the input value is
		 * returned unmodified.
		 */
		double GetPeriodicValue(double Location) const override
		{
			if(!_periodic)
				return Location;

			int n = (int)(Location/(2.0*M_PI));
			double PeriodicLocation = Location-2.0*n*M_PI;

			PeriodicLocation = Location - n*M_PI;
			if(PeriodicLocation < -M_PI)
				PeriodicLocation += 2.0*M_PI;
			else if (Location > M_PI)
				PeriodicLocation -= 2.0*M_PI;

			return PeriodicLocation;
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
		 * \param Location Input value
		 * \return Wrapped or unwrapped difference depending on whether periodic
		 *         boundaries are used.
		 *
		 * If periodic boundaries are used, this function calculates the
		 * difference and wraps the result into the range (-pi, pi). Otherwise,
		 * the simple difference is returned.
		 */
		double GetDifference(const double Location) const override
		{
			double PeriodicDiff = _val - Location;

			if(!_periodic)
				return PeriodicDiff;

			PeriodicDiff = GetPeriodicValue(PeriodicDiff);

			if(PeriodicDiff > M_PI)
				PeriodicDiff -= 2.0*M_PI;
			else if(PeriodicDiff < -M_PI)
				PeriodicDiff += 2.0*M_PI;

			return PeriodicDiff;
		}

		//! Serialize this CV for restart purposes.
		/*!
		 * \param json JSON value
		 */
		virtual void Serialize(Json::Value& json) const override
		{
			json["type"] = "Torsional";
			json["periodic"] = _periodic;
			json["atom ids"][0] = _atomid1;
			json["atom ids"][1] = _atomid2;
			json["atom ids"][2] = _atomid3;
			json["atom ids"][3] = _atomid4;
			for(size_t i = 0; i < _bounds.size(); ++i)
				json["bounds"].append(_bounds[i]);

		}
	};
}

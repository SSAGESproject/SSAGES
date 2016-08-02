#pragma once 

#include "CollectiveVariable.h"

#include <array>
#include <math.h>
#include "../Utility/VectorProducts.h"

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
			double ix = 0;
			double iy = 0;
			double iz = 0;
			double jx = 0;
			double jy = 0;
			double jz = 0;
			double kx = 0;
			double ky = 0;
			double kz = 0;
			double lx = 0;
			double ly = 0;
			double lz = 0;
			
			// Loop through atom positions
			for(size_t i = 0; i < pos.size(); ++i)
			{
				_grad[i][0] = 0;
				_grad[i][1] = 0;
				_grad[i][2] = 0;
				// If we are at the atom ID of interest, grab coordinates
				if(ids[i] == _atomid1)
				{
					//coordinates for atom i
					ix = pos[i][0];
					iy = pos[i][1];
					iz = pos[i][2];
					iindex = i;
				}
				if(ids[i] == _atomid2)
				{
					//coordinates for atom j
					jx = pos[i][0];
					jy = pos[i][1];
					jz = pos[i][2];
					jindex = i;
				}
				if(ids[i] == _atomid3)
				{
					//coordinates for atom k
					kx = pos[i][0];
					ky = pos[i][1];
					kz = pos[i][2];
					kindex = i;
				}
				if(ids[i] == _atomid4)
				{
					//coordinates for atom l
					lx = pos[i][0];
					ly = pos[i][1];
					lz = pos[i][2];
					lindex = i;
				}
			}

			if(iindex < 0 || jindex < 0 || kindex < 0 || lindex <0)
			{
				std::cout<<"Out of bounds index, could not locate an ID"<<std::endl;
				exit(0);
			}

			//Calculate pertinent vectors
			Vector3 F{
				ix - jx,
				iy - jy,
				iz - jz};
			Vector3 G{
				jx - kx,
				jy - ky,
				jz - kz};
			Vector3 H{
				lx - kx,
				ly - ky, 
				lz - kz};

			//Vector3 rij{{
			//	ix - jx,
			//	iy - jy,
			//	iz - jz}};
			//Vector3 rkj{{
			//	kx - jx,
			//	ky - jy,
			//	kz - jz}};
			//Vector3 rkl{{
			//	lx - kx,
			//	ly - ky, 
			//	lz - kz}};

			//Calculate dihedral angle

			//Vector3 rim, rln;
			//double rkj2 = norm(rkj)*norm(rkj);
			//for(size_t i = 0; i<3; i++)
		//	{
		//		rim[i] = rij[i] - DotProduct(rij,rkj)/(rkj2)*rkj[i];
		//		rln[i] = DotProduct(rkl,rkj)/(rkj2)*rkj[i] - rkl[i];
		//	}
//
//			double normrim, normrln;
//
//			normrim = norm(rim);
//			normrln = norm(rln);
//
//			auto normkj = norm(rkj);
//			Vector3 rijrkjprod;
//			for(size_t i = 0; i < rijrkjprod.size();i++)
//				rijrkjprod[i] = rij[i]*normkj;
//			auto rkjrklcross = CrossProduct(rkj, rkl);
//			auto y = DotProduct(rijrkjprod, rkjrklcross);
//			auto rijrkjcross = CrossProduct(rij, rkj);
//			auto x = DotProduct(rijrkjcross, rkjrklcross);
//			_val = atan2(y, x);

			Vector3 A = CrossProduct(F, G);
			Vector3 B = CrossProduct(H, G);

			auto y1 = CrossProduct(B, A);
			auto y = DotProduct(y1, G)/norm(G);
			auto x = DotProduct(A, B);

			_val = atan2(y, x);


			double Zed = DotProduct(F, G)/(norm2(A)*norm(G));
			double Ned = DotProduct(H, G)/(norm2(B)*norm(G));

			for(size_t i = 0; i<3; i++)
			{
				_grad[iindex][i] = -norm(G)*A[i]/norm2(A);
				_grad[lindex][i] = norm(G)*B[i]/norm2(B);
				_grad[jindex][i] = Zed*A[i] - Ned*B[i] - _grad[iindex][i];
				_grad[kindex][i] = Ned*B[i] - Zed*A[i] - _grad[lindex][i];
			}

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

			double pi = 3.14159;
			int n = (int)(Location/(2.0*pi));
			double PeriodicLocation = Location-2.0*n*pi;

			PeriodicLocation = Location - n*pi;
			if(PeriodicLocation < -pi)
				PeriodicLocation += 2.0*pi;
			else if (Location > pi)
				PeriodicLocation -= 2.0*pi;

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
			double pi = 3.14159;
			double PeriodicDiff = _val - Location;

			if(!_periodic)
				return PeriodicDiff;

			PeriodicDiff = GetPeriodicValue(PeriodicDiff);

			if(PeriodicDiff > pi)
				PeriodicDiff -= 2.0*pi;
			else if(PeriodicDiff < -pi)
				PeriodicDiff += 2.0*pi;

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

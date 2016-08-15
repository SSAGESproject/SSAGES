#pragma once 

#include "CollectiveVariable.h"

#include <array>
#include <cmath>

namespace SSAGES
{
	// Collective variable to calculate radius of gyration.

	class AngleCV : public CollectiveVariable
	{
	private:
		// ID of first atom
		int _atomid1;

		// ID of second atom 
		int _atomid2;

		// ID of third atom
		int _atomid3;

		// Current value of the CV.
		double _val;

		// Gradient of the CV, dRg/dxi.
		std::vector<Vector3> _grad;

		// Bounds on CV.
		std::array<double, 2> _bounds;

	public:

		//! Constructor.
		/*!
		 * \param atomid1 ID of the first atom defining the angle.
		 * \param atomid2 ID of the second atom defining the angle.
		 * \param atomid3 ID of the third atom defining the angle.
		 * \param periodic If \c True consider periodic boundary conditions.
		 *
		 * Construct an dihedral CV.
		 *
		 * \todo Bounds needs to be an input and periodic boundary conditions
		 */
		AngleCV(int atomid1, int atomid2, int atomid3) :
		_atomid1(atomid1),_atomid2(atomid2), _atomid3(atomid3), _val(0), _grad(0), _bounds{{0,M_PI}}
		{
		}

		//! Initialize necessary variables.
		/*!
		 * \param snapshot Current simulation snapshot.
		 */
		void Initialize(const Snapshot& snapshot) override
		{
			// Initialize gradient
			auto n = snapshot.GetPositions().size();
			_grad.resize(n);

		}

		//! Evaluate the CV.
		/*!
		 * \param snapshot Current simulation snapshot.
		 */
		void Evaluate(const Snapshot& snapshot) override
		{
			const auto& pos = snapshot.GetPositions();
			const auto& ids = snapshot.GetAtomIDs();

			Vector3 atomi, atomj, atomk;
			int iindex = -1, jindex = -1, kindex = -1;

			// Loop through atom positions
			for( size_t i = 0; i < pos.size(); ++i)
			{
				_grad[i].setZero();
				
				if( ids[i] == _atomid1)
				{
					atomi = pos[i];
					iindex=i;
				}
				else if( ids[i] == _atomid2)
				{
					atomj = pos[i];
					jindex=i;
				}
				else if(ids[i] == _atomid3)
				{
					atomk = pos[i];
					kindex=i;
				}
			}

			if(iindex < 0 || jindex < 0 || kindex < 0 )
			{
				std::cout<<"Out of bounds index, could not locate an ID"<<std::endl;
				exit(0);
			}

			// Two vectors
			Vector3 rij = snapshot.ApplyMinimumImage(atomi - atomj);
			Vector3 rkj = snapshot.ApplyMinimumImage(atomk - atomj);

			auto dotP = rij.dot(rkj); //DotProduct(rij, rkj);
			auto nrij = rij.norm(); //norm(rij);
			auto nrkj = rkj.norm(); //norm(rkj);

			_val = acos(dotP/(nrij*nrkj));

			// Calculate gradients
			double prefactor = -1.0/(sqrt(1 - dotP/(nrij*nrkj))*nrij*nrkj);

			_grad[iindex] = prefactor*(rkj - dotP*rij/(nrij*nrij));	
			_grad[kindex] = prefactor*(rij - dotP*rkj/(nrkj*nrkj));
			_grad[jindex] = -_grad[iindex] - _grad[kindex];
		}

		//! Return the value of the CV.
		/*!
		 * \return Current value of the CV.
		 */
		double GetValue() const override 
		{ 
			return _val; 
		}

		//! Return value taking periodic boundary conditions into account.
		/*!
		 * \param Location Input value.
		 * \return Original value.
		 *
		 * The atom position CV does not take periodic boundary conditions into
		 * account. Thus, this function will always return the unmodified input
		 * value.
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
		 * \return Values of the lower and upper boundaries of this CV.
		 */
		const std::array<double, 2>& GetBoundaries() const override
		{
			return _bounds;
		}

		//! Get difference taking periodic boundary conditions into account.
		/*!
		 * \param Location Input value.
		 * \return Simple difference between current value of CV and Location.
		 *
		 * The atom position CV does not take periodic boundary conditions into
		 * account. Thus, this function returns the simple difference between
		 * the current value of the CV and the input value.
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
			json["type"] = "Angle";
			json["atom ids"][0] = _atomid1;
			json["atom ids"][1] = _atomid2;
			json["atom ids"][2] = _atomid3;
			for(size_t i = 0; i < _bounds.size(); ++i)
				json["bounds"].append(_bounds[i]);

		}

	};
}

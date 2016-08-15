#pragma once 

#include "CollectiveVariable.h"

#include <array>
#include <math.h>

namespace SSAGES
{
	//! Collective variable on an atom position.
	/*!
	 * This CV will return the distance of an atom from a particular point in
	 * (1,2,3)-dimensional space.
	 *
	 * \ingroup CVs
	 */
	class AtomPositionCV : public CollectiveVariable
	{
	private:
		//! ID of atom of interest.
		int _atomid; 

		//! Point in space.
		Vector3 _position;

		//! Current value of the CV.
		double _val;

		bool _fixx; //!< If \c False, constrain x dimension of CV.
		bool _fixy; //!< If \c False, constrain y dimension of CV.
		bool _fixz; //!< If \c False, constrain z dimension of CV.

		//! Gradient of the CV, dr/dxi.
		std::vector<Vector3> _grad;

		//! Bounds on CV.
		std::array<double, 2> _bounds;

	public:
		//! Constructor
		/*!
		 * \param atomid ID of the atom of interest.
		 * \param position Point in (1,2,3)D space for the distance calculation.
		 * \param fixx If \c False, constrain x dimension.
		 * \param fixy If \c False, constrain y dimension.
		 * \param fixz If \c False, constrain z dimension.
		 *
		 * Construct an atom position CV. If a dimension is constrained, this
		 * dimension does not contribute to the distance calculation. For
		 * example, if the y- and z-dimension are constrained, the CV calculates
		 * only the distance in x-direction.
		 *
		 * \todo Bounds needs to be an input.
		 */
		AtomPositionCV(int atomid, const Vector3& position, bool fixx, bool fixy, bool fixz) : 
		_atomid(atomid), _position(position), _val(0), 
		_fixx(fixx), _fixy(fixy), _fixz(fixz), _grad(0), _bounds{{0,0}}
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
			
			// Loop through atom positions
			for(size_t i = 0; i < pos.size(); ++i)
			{
				// If we are at the atom ID of interest.
				if(ids[i] == _atomid)
				{
					// Compute distance.
					Vector3 dx = snapshot.ApplyMinimumImage(pos[i] - _position);

					// Set to 0 dimensions we aren't interested in.
					if(!_fixx) dx[0] = 0;
					if(!_fixy) dx[1] = 0;
					if(!_fixz) dx[2] = 0;

					// Compute norm and gradient.
					auto r = dx.norm();
					_val = r;
					if (r == 0)
						_grad[i].setZero();
					else 
						_grad[i] = dx/r;
				}
				else
				{
					_grad[i].setZero();
				}
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

		}
	};
}

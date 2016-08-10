#pragma once 

#include "CollectiveVariable.h"

#include <array>
#include <math.h>

namespace SSAGES
{
	//! Collective variable on an atom coordinate.
	/*!
	 * This will return the value of either the x, y, or z coordinate, depending
	 * on the user specification for a defined atom.
	 *
	 * \ingroup CVs
	 */
	class AtomCoordinateCV : public CollectiveVariable
	{
	private:
		//! ID of atom of interest.
		int _atomid; 

		//! Index of dimension. 0 -> x, 1 -> y, 2 -> z.
		int _index;

		//! Current value of CV.
		double _val;

		//! Gradient vector dOP/dxi.
		std::vector<Vector3> _grad;

		//! Bounds on CV.
		std::array<double, 2> _bounds;

	public:
		//! Constructor.
		/*!
		 * \param atomid Index of the atom of interest.
		 * \param index Index of dimension.
		 *
		 * Construct an atom coordinate CV. The atomid specifies the ID of the
		 * atom of interest, and index specifies the dimension to report with
		 * 0 -> x, 1 -> y, 2 -> z.
		 *
		 * \todo Bounds needs to be an input.
		 */
		AtomCoordinateCV(int atomid, int index) : 
		_atomid(atomid), _index(index), _val(0), _grad(0), _bounds{{0,0}}
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
			const auto& ids = snapshot.GetAtomIDs();
			_grad.resize(n);

			if(std::find(ids.begin(), ids.end(), _atomid) == ids.end())
			    throw BuildException({"AtomCoordinateCV: Cannot find specified ID."});

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

			// Loop through atom positions.
			for(size_t i = 0; i < ids.size(); ++i)
			{
				_grad[i].setZero();
				// We are at the ID of interest.
				if(ids[i] == _atomid)
				{
					// Set to unity the dimension of interest.
					switch(_index)
					{
						case 0:
							_val = pos[i][0];
							_grad[i][0] = 1.0;
							break;
						case 1:
							_val = pos[i][1];
							_grad[i][1] = 1.0;
							break;
						case 2:
							_val = pos[i][2];
							_grad[i][2] = 1.0;
							break;
					}
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
		 * \return Input value
		 *
		 * The AtomCoordinate CV does not consider periodic boundary
		 * conditions. Thus, this function always returns the input value.
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
		 * \return Values of the lower and upper boundary.
		 */
		const std::array<double, 2>& GetBoundaries() const override
		{
			return _bounds;
		}

		//! Return difference considering periodic boundaries.
		/*!
		 * \return Direct difference.
		 *
		 * As the AtomCoordinate CV does not consider periodic boundary
		 * conditions, the difference between the current value of the CV and
		 * another value is always the direct difference.
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
			json["type"] = "AtomCoordinate";			
			switch(_index)
					{
						case 0:
							json["dimension"] = "x";
							break;
						case 1:
							json["dimension"] = "y";
							break;
						case 2:
							json["dimension"] = "z";
							break;
					}
			json["atom id"] = _atomid;
			for(size_t i = 0; i < _bounds.size(); ++i)
				json["bounds"].append(_bounds[i]);

		}
	};
}

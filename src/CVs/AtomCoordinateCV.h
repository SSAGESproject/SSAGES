#pragma once 

#include "CollectiveVariable.h"

#include <array>
#include <math.h>

namespace SSAGES
{
	// Collective variable on an atom coordinate. This will
	// return the value of either the x, y, or z coordinate
	// depending on the user specification for a defined atom.
	class AtomCoordinateCV : public CollectiveVariable
	{
	private:
		// ID of atom of interest.
		int _atomid; 

		// Index of dimension. 0 -> x, 1 -> y, 2 -> z.
		int _index;

		// Current value of CV.
		double _val;

		// Gradient vector dOP/dxi.
		std::vector<Vector3> _grad;

		// Bounds on CV.
		std::array<double, 2> _bounds;

	public:
		// Construct an atom coordinate CV. The atomid specifies the 
		// ID of the atom of interest, and index specifies the dimension 
		// to report with 0 -> x, 1 -> y, 2 -> z. 
		// TODO: bounds needs to be an input.
		AtomCoordinateCV(int atomid, int index) : 
		_atomid(atomid), _index(index), _val(0), _grad(0), _bounds{{0,0}}
		{
		}

		// Initialize necessary variables.
		void Initialize(const Snapshot& snapshot) override
		{
			// Initialize gradient. 
			auto n = snapshot.GetPositions().size();
			_grad.resize(n);
		}

		// Evaluate the CV.
		void Evaluate(const Snapshot& snapshot) override
		{
			// Gradient and value. 
			const auto& pos = snapshot.GetPositions(); 
			const auto& ids = snapshot.GetAtomIDs();

			// Loop through atom positions.
			for(size_t i = 0; i < pos.size(); ++i)
			{
				// We are at the ID of interest.
				if(ids[i] == _atomid)
				{
					// We set the gradient to zero, and only 
					// set to unity the dimension of interest.
					_grad[i][0] = 0;
					_grad[i][1] = 0;
					_grad[i][2] = 0;
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
				else
				{
					_grad[i][0] = 0;
					_grad[i][1] = 0;
					_grad[i][2] = 0;
				}
			}
		}

		// Return the value of the CV.
		double GetValue() const override 
		{ 
			return _val; 
		}

		double GetPeriodicValue(double Location) const override
		{
			return Location;
		}

		// Return the gradient of the CV.
		const std::vector<Vector3>& GetGradient() const override
		{
			return _grad;
		}

		// Return the boundaries of the CV.
		const std::array<double, 2>& GetBoundaries() const override
		{
			return _bounds;
		}

		double GetDifference(const double Location) const override
		{
			return _val - Location;
		}
	};
}

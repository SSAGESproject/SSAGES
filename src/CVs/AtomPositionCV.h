#pragma once 

#include "CollectiveVariable.h"

#include <array>
#include <math.h>

namespace SSAGES
{
	// Collective variable on an atom position. This will 
	// return the distance of an atom from a particular point 
	// in (1,2,3)-dimensional space.
	class AtomPositionCV : public CollectiveVariable
	{
	private:
		// ID of atom of interest.
		int _atomid; 

		// Point in space.
		Vector3 _position;

		// Current value of the CV.
		double _val;

		// Constraints in x,y,z dimensions.
		bool _fixx, _fixy, _fixz;

		// Gradient of the CV, dr/dxi.
		std::vector<Vector3> _grad;

		// Bounds on CV.
		std::array<double, 2> _bounds;

		// Helper function to compute the norm of a vector.
		double norm(const Vector3& v)
		{
			return sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
		}

	public:
		// Construct an atom position CV. The atomid specifies 
		// the ID of the atom of interest and position is the 
		// point in (1,2,3)D space for the distance calculation.
		// The fix booleans specify if that dimension is to be 
		// constrained.
		// TODO: bounds needs to be an input.
		AtomPositionCV(int atomid, const Vector3& position, bool fixx, bool fixy, bool fixz) : 
		_atomid(atomid), _position(position), _val(0), 
		_fixx(fixx), _fixy(fixy), _fixz(fixz), _grad(0), _bounds{{0,0}}
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
			
			// Loop through atom positions
			for(size_t i = 0; i < pos.size(); ++i)
			{
				// If we are at the atom ID of interest.
				if(ids[i] == _atomid)
				{
					// Compute distance.
					Vector3 dx{{
						pos[i][0] - _position[0], 
						pos[i][1] - _position[1], 
						pos[i][2] - _position[2]}};

					// Set to 0 dimensions we aren't interested in.
					if(!_fixx) dx[0] = 0;
					if(!_fixy) dx[1] = 0;
					if(!_fixz) dx[2] = 0;

					// Compute norm and gradient.
					auto r = norm(dx);
					_grad[i][0] = dx[0]/r;
					_grad[i][1] = dx[1]/r;
					_grad[i][2] = dx[2]/r;
					_val = r;
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
	};
}
#pragma once 

#include "CollectiveVariable.h"

#include <array>
#include <math.h>


namespace SSAGES
{
	// Collective variable on an atom position. This will 
	// return the distance of an atom from a particular point 
	// in (1,2,3)-dimensional space.
	class AtomSeparationCV : public CollectiveVariable
	{
	private:
		// ID of atom of interest 1.
		int _atomid1;
		
		// ID of atom of interest 1.
		int _atomid2;

		// Current value of the CV.
		double _val;

		// Simulation box size.
		std::vector<double> _boxsize;

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
		AtomSeparationCV(int atomid1, int atomid2, std::vector<double> boxsize) : 
		_atomid1(atomid1),_atomid2(atomid2), _boxsize(boxsize), _val(0), _grad(0), _bounds{{0,0}}
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
			
			// Some temp variables.
			//double posx1 = 0,posx2 = 0,posy1 = 0,posy2 = 0,posz1 = 0,posz2 = 0;
			Vector3 pos1{{0, 0, 0}};
			Vector3 pos2{{0, 0, 0}};
			
			size_t index1 = 0, index2 = 0;
		
			
			// Loop through atom positions
			for(size_t i = 0; i < pos.size(); ++i)
			{
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
				else
				{
					_grad[i][0] = 0;
					_grad[i][1] = 0;
					_grad[i][2] = 0;
				}
			}

			Vector3 del{{pos1[0]-pos2[0],pos1[1]-pos2[1],pos1[2]-pos2[2]}};
			
			for(size_t i = 0; i<del.size(); ++i)
				{
				while(del[i] < -_boxsize[i]/2){
					del[i] += _boxsize[i];
					}
				while(del[i] > _boxsize[i]/2){
					del[i] -= _boxsize[i];
					}
				}					

			// Compute norm and gradient.
			auto r = norm(del);
			for(size_t i = 0; i<del.size(); ++i)
				{
				_grad[index1][i] = del[i]/r;
				_grad[index2][i] = -del[i]/r;
				}
			_val = r;
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

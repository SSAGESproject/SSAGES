#pragma once 

#include "CollectiveVariable.h"

#include <array>
#include <math.h>

namespace SSAGES
{
	class AtomCoordinateCV : public CollectiveVariable
	{
	private:
		int _atomid; 
		int _index;
		double _val;
		std::vector<Vector3> _grad;
		std::array<double, 2> _bounds;

		double norm(const Vector3& v)
		{
			return sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
		}

	public:
		AtomCoordinateCV(int atomid, int index) : 
		_atomid(atomid), _index(index), _val(0), _grad(0), _bounds{{0,0}}
		{
		}

		void Initialize(const Snapshot& snapshot) override
		{
			// Initialize gradient. 
			auto n = snapshot.GetPositions().size();
			_grad.resize(n);
		}

		void Evaluate(const Snapshot& snapshot) override
		{
			// Gradient and value. 
			auto& pos = snapshot.GetPositions(); 
			auto& ids = snapshot.GetAtomIDs();
			for(size_t i = 0; i < pos.size(); ++i)
			{
				if(ids[i] == _atomid)
				{
					//grab positions, decide which one of the three directions matters, and move onward.
					double dx = 0;
					_grad[i][0] = 0;
					_grad[i][1] = 0;
					_grad[i][2] = 0;
					switch(_index)
					{
						case 0:
							dx = pos[i][0];
							_grad[i][0] = 1.0;
							break;
						case 1:
							dx = pos[i][1];
							_grad[i][1] = 1.0;
							break;
						case 2:
							dx = pos[i][2];
							_grad[i][2] = 1.0;
							break;
					}
					_val = dx;
				}
				else
				{
					_grad[i][0] = 0;
					_grad[i][1] = 0;
					_grad[i][2] = 0;
				}
			}
		}

		double GetValue() const override 
		{ 
			return _val; 
		}

		const std::vector<Vector3>& GetGradient() const override
		{
			return _grad;
		}

		const std::array<double, 2>& GetBoundaries() const override
		{
			return _bounds;
		}
	};
}

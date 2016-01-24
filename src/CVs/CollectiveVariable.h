#pragma once
#include "../Snapshot.h"

#include <vector>

namespace SSAGES
{
	class CollectiveVariable
	{
	public:
		virtual void Initialize(const Snapshot&) {}
		virtual void Evaluate(const Snapshot& snapshot) = 0;
		virtual double GetValue() const = 0;
		virtual const std::vector<Vector3>& GetGradient() const = 0;
		virtual const std::array<double, 2>& GetBoundaries() const = 0;
		//virtual const std::vector<double>& GetAtomIDs() const = 0;
	};

	// Definitions.
	using CVList = std::vector<CollectiveVariable*>;
}
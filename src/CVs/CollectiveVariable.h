#pragma once
#include "../Snapshot.h"

#include <vector>

// Forward declare.
namespace Json {
	class Value;
}

namespace SSAGES
{
	// Interface for a collective variable.
	class CollectiveVariable
	{
	public:
		// Initialization of the CV. This is an optional 
		// method and is called during the pre-simulation phase
		// of the hook. It is typically used to allocate/reserve
		// memory. 
		virtual void Initialize(const Snapshot&) {}

		// Evaluation of the CV. This is called by the Hook in the 
		// post-integration phase every iteration. The CV should compute
		// its value and gradient, storing them in a local private variable.
		virtual void Evaluate(const Snapshot& snapshot) = 0;

		// Returns the current value of the CV which has been computed 
		// via the call to Evaluate. This is const, and is typically
		// used by Methods.
		virtual double GetValue() const = 0;

		// Takes location and applies periodic boundaries of the CV on it
		// and returns a correct value. Example would be torsional angle
		// which has bounds at pi and -pi. If location = 2pi,
		// GetPeriodicValue(location) would return 0.
		virtual double GetPeriodicValue(double Location) const = 0;

		// Returns the current value of the CV gradient. This should be 
		// an n length vector, where n is the number of atoms in the 
		// snapshot. Each element in the vector is the derivative of the 
		// CV w.r.t. the atom's coordinate (dCV/dxi).
		virtual const std::vector<Vector3>& GetGradient() const = 0;

		// Returns the boundaries of the CV. These represent the bounds 
		// within which the CV is expected to be constrained. There is 
		// no requirement on the method to respect the values returned here.
		virtual const std::array<double, 2>& GetBoundaries() const = 0;

		// Returns the difference betwen the current cv value and Location:
		// (_value - Location) respecting periodic boundary conditions of the CV,
		// if the CV has periodic boundary conditions. 
		// For example Torsional angle has boundaries at
		// pi and -pi, in which the difference beteen the angles is 0 not 2pi
		virtual double GetDifference(const double Location) const = 0;

		// Builds a cv from a JSON node. Returns a pointer to the built cv.
		// If return value is nullptr, 
		// then an unknown error occurred. It will throw a BuildException on failure. 
		// Object lifetime is the caller's responsibility. 
		static CollectiveVariable* BuildCV(const Json::Value& json);

		// Overloaded function allowing JSON path specification.
		static CollectiveVariable* BuildCV(const Json::Value& json, 
							   const std::string& path);
	};

	// Definitions.
	using CVList = std::vector<CollectiveVariable*>;
}
#pragma once 

#include "Method.h"
#include "../CVs/CollectiveVariable.h"

namespace SSAGES
{
	class Umbrella : public Method
	{
	private:
		std::vector<double> _kspring;
		std::vector<double> _centers;

	public:
		Umbrella(const std::vector<double>& kspring,
				 const std::vector<double>& centers,
				 unsigned int frequency) : 
		Method(frequency), _kspring(kspring), _centers(centers)
		{}

		// Pre-simulation hook.
		void PreSimulation(Snapshot* snapshot, const CVList& cvs) override;

		// Post-integration hook.
		void PostIntegration(Snapshot* snapshot, const CVList& cvs) override;

		// Post-simulation hook.
		void PostSimulation(Snapshot* snapshot, const CVList& cvs) override;
	};
}
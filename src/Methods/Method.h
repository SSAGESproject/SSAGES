#pragma once

#include "../EventListener.h"

namespace SSAGES
{
	// Abstract base class for method implementations.
	class Method : public EventListener
	{
	public:
		Method(unsigned int frequency) : 
		EventListener(frequency) {}

		// Method call prior to simulation initiation.
		virtual void PreSimulation(Snapshot* snapshot, const CVList& cvs) override = 0;

		// Method call post integration.
		virtual void PostIntegration(Snapshot* snapshot, const CVList& cvs) override = 0;

		// Method call post simulation.
		virtual void PostSimulation(Snapshot* snapshot, const CVList& cvs) override = 0;

		~Method() {}
	};
}
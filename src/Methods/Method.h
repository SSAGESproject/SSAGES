#pragma once

#include "../EventListener.h"

namespace SSAGES
{
	// Interface for Method implementations.
	class Method : public EventListener
	{
	public:
		// Frequency of sampling must be specified by all methods.
		Method(unsigned int frequency) : 
		EventListener(frequency) {}

		// Method call prior to simulation initiation.
		virtual void PreSimulation(Snapshot* snapshot, const CVList& cvs) override = 0;

		// Method call post integration.
		virtual void PostIntegration(Snapshot* snapshot, const CVList& cvs) override = 0;

		// Method call post simulation.
		virtual void PostSimulation(Snapshot* snapshot, const CVList& cvs) override = 0;

		virtual ~Method() {}
	};
}
# pragma once 

#include "Snapshot.h"

namespace SSAGES
{
	// Base abstract class for listening in to events fired by "Hook".
	class EventListener
	{
	private:
		unsigned int _frequency;
	
	public:
		EventListener(unsigned int frequency) : 
		_frequency(frequency)
		{

		}

		unsigned int GetFrequency() const { return _frequency; }

		// Method call prior to simulation initiation.
		virtual void PreSimulation(Snapshot* snapshot) = 0;

		// Method call post integration.
		virtual void PostIntegration(Snapshot* snapshot) = 0;

		// Method call post simulation.
		virtual void PostSimulation(Snapshot* snapshot) = 0;
	};
}
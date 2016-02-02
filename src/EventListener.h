# pragma once 

#include "Snapshot.h"
#include "CVs/CollectiveVariable.h"

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

		// Get frequency of event listener.
		unsigned int GetFrequency() const { return _frequency; }

		// Method call prior to simulation initiation.
		virtual void PreSimulation(Snapshot* snapshot, const CVList& cvs) = 0;

		// Method call post integration.
		virtual void PostIntegration(Snapshot* snapshot, const CVList& cvs) = 0;

		// Method call post simulation.
		virtual void PostSimulation(Snapshot* snapshot, const CVList& cvs) = 0;
	};
}
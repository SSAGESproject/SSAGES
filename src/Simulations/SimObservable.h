#pragma once

#include "../Observers/Visitable.h"
#include "SimEvent.h"
#include "SimObserver.h"
#include <list>

namespace SSAGES
{
	// Abstract class for observables in a simulation.
	class SimObservable : public Visitable
	{
	private:
		std::vector<SimObserver*> _observers;

	protected:
		void SerializeObservers(Json::Value& json) const
		{
			for(int i = 0; i < (int)_observers.size(); ++i)
				_observers[i]->Serialize(json[i]);
		}

	public:

		// Add simulation observer.
		virtual void AddObserver(SimObserver* observer);

		// Remove simulation observer.
		virtual void RemoveObserver(SimObserver* observer);

		// Notify registered observers of a change.
		void NotifyObservers(const SimEvent& event);
	};
}

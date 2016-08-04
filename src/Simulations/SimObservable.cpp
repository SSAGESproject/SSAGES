#include "SimObservable.h"
#include <algorithm>

namespace SSAGES
{
	// Add simulation observer.
	void SimObservable::AddObserver(SimObserver* observer)
	{
		_observers.push_back(observer);
	}

	// Remove simulation observer.
	void SimObservable::RemoveObserver(SimObserver* observer)
	{
		_observers.erase(
			std::remove(_observers.begin(), _observers.end(), observer), 
			_observers.end()
		);
	}

	// Notify all observers of a simulation event.
	void SimObservable::NotifyObservers(const SimEvent& event)
	{
		for(auto& observer : _observers)
			observer->Update(event);
	}
}

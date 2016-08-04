#pragma once

#include "../Observers/Visitable.h"
#include "SimEvent.h"
#include "SimObserver.h"
#include <list>

namespace SSAGES
{
	//! Base class for observable objects in a simulation.
	/*
	 * Abstract class for observables in a simulation.
	 *
	 * \ingroup Core
	 */
	class SimObservable : public Visitable
	{
	private:
		//! List of observers observing this object.
		std::vector<SimObserver*> _observers;

	protected:
		//! Serialize observers
		/*
		 * \param json JSON value to store observer information to.
		 */
		void SerializeObservers(Json::Value& json) const
		{
			for(int i = 0; i < (int)_observers.size(); ++i)
				_observers[i]->Serialize(json[i]);
		}

	public:

		//! Add simulation observer.
		/*!
		 * \param observer Observer to add.
		 */
		virtual void AddObserver(SimObserver* observer);

		//! Remove simulation observer.
		/*!
		 * \param observer Pointer to observer to be removed.
		 */
		virtual void RemoveObserver(SimObserver* observer);

		//! Notify registered observers of a change.
		/*!
		 * \param event Simulation event that has happened.
		 */
		void NotifyObservers(const SimEvent& event);

		//! Destructor.
		virtual ~SimObservable() {}
	};
}

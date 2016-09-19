/**
 * This file is part of
 * SSAGES - Suite for Advanced Generalized Ensemble Simulations
 *
 * Copyright 2016 Ben Sikora <bsikora906@gmail.com>
 *
 * SSAGES is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * SSAGES is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with SSAGES.  If not, see <http://www.gnu.org/licenses/>.
 */
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

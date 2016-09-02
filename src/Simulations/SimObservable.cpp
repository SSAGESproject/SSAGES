/**
 * This file has been adapted from
 * SAPHRON - Statistical Applied PHysics through Random On-the-fly Numerics
 * https://github.com/hsidky/SAPHRON
 *
 * Copyright 2016 Hythem Sidky
 *
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

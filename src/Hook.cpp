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

#include "Hook.h"
#include "Drivers/Driver.h"
#include <algorithm>

namespace SSAGES
{
	//! Sets the active snapshot.
	void Hook::SetSnapshot(Snapshot* snapshot)
	{
		_snapshot = snapshot;
	}

	//! Sets the active Driver
	void Hook::SetMDDriver(Driver* MDDriver)
	{
		_MDDriver = MDDriver;
	}

	//! Add a listener to the hook.
	/*!
	 * \param listener Pointer to the EventListener to be added to the Hook.
	 *
	 * Does nothing if the listener is already added.
	 */
	void Hook::AddListener(EventListener* listener)
	{
		if(std::find(_listeners.begin(), _listeners.end(), listener) == _listeners.end())
			_listeners.push_back(listener);
	}

	//! Add a CollectiveVariable to the hook.
	/*!
	 * \param cv CollectiveVariable to be added to the hook.
	 *
	 * Does nothing if the CollectiveVariable is already added.
	 */
	void Hook::AddCV(CollectiveVariable* cv)
	{
		if(std::find(_cvs.begin(), _cvs.end(), cv) == _cvs.end())
			_cvs.push_back(cv);
	}

	void Hook::NotifyObservers()
	{
		auto& comm = _snapshot->GetCommunicator();
		if(comm.rank() == 0)
			_MDDriver->NotifyObservers(SimEvent(_MDDriver, _snapshot->GetIteration()));
	}
}
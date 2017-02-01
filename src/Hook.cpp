/**
 * This file is part of
 * SSAGES - Suite for Advanced Generalized Ensemble Simulations
 *
 * Copyright 2016 Ben Sikora <bsikora906@gmail.com>
 *                Hythem Sidky <hsidky@nd.edu>
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
#include "CVs/CollectiveVariable.h"
#include <algorithm>

namespace SSAGES
{
	void Hook::PreSimulationHook()
	{
		_snapshot->Changed(false);
	
		// Initialize/evaluate CVs.
		for(auto& cv : _cvs)
		{
			cv->Initialize(*_snapshot);
			cv->Evaluate(*_snapshot);
		}

		// Call presimulation method on listeners. 
		for(auto& listener : _listeners)
			listener->PreSimulation(_snapshot, _cvs);

		// Sync snapshot to engine.
		if(_snapshot->HasChanged())
			SyncToEngine();

		_snapshot->Changed(false);		
	}

	void Hook::PostIntegrationHook()
	{
		_snapshot->Changed(false);

		for(auto& cv : _cvs)
			cv->Evaluate(*_snapshot);

		for(auto& listener : _listeners)
			if(_snapshot->GetIteration() % listener->GetFrequency() == 0)
				listener->PostIntegration(_snapshot, _cvs);

		if(_snapshot->HasChanged())
			SyncToEngine();

		_snapshot->Changed(false);		
	}

	void Hook::PostSimulationHook()
	{
		_snapshot->Changed(false);

		for(auto& cv : _cvs)
			cv->Evaluate(*_snapshot);
		
		for(auto& listener : _listeners)
			listener->PostSimulation(_snapshot, _cvs);

		if(_snapshot->HasChanged())
			SyncToEngine();

		_snapshot->Changed(false);		
	}

	//! Sets the active snapshot.
	void Hook::SetSnapshot(Snapshot* snapshot)
	{
		snapshot_ = snapshot;
	}

	//! Sets the active Driver
	void Hook::SetMDDriver(Driver* MDDriver)
	{
		MDDriver_ = MDDriver;
	}

	//! Add a listener to the hook.
	/*!
	 * \param listener Pointer to the EventListener to be added to the Hook.
	 *
	 * Does nothing if the listener is already added.
	 */
	void Hook::AddListener(EventListener* listener)
	{
		if(std::find(listeners_.begin(), listeners_.end(), listener) == listeners_.end())
			listeners_.push_back(listener);
	}

	//! Add a CollectiveVariable to the hook.
	/*!
	 * \param cv CollectiveVariable to be added to the hook.
	 *
	 * Does nothing if the CollectiveVariable is already added.
	 */
	void Hook::AddCV(CollectiveVariable* cv)
	{
		if(std::find(cvs_.begin(), cvs_.end(), cv) == cvs_.end())
			cvs_.push_back(cv);
	}

	void Hook::NotifyObservers()
	{
		auto& comm = snapshot_->GetCommunicator();
		if(comm.rank() == 0)
			MDDriver_->NotifyObservers(SimEvent(MDDriver_, snapshot_->GetIteration()));
	}
}
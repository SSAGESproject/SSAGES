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
		snapshot_->Changed(false);
	
		// Initialize/evaluate CVs.
		for(auto& cv : cvmanager_->GetCVs({}))
		{
			cv->Initialize(*snapshot_);
			cv->Evaluate(*snapshot_);
		}

		// Call presimulation method on listeners. 
		for(auto& listener : listeners_)
			listener->PreSimulation(snapshot_, cvmanager_);

		// Sync snapshot to engine.
		if(snapshot_->HasChanged())
			SyncToEngine();

		snapshot_->Changed(false);		
	}

	void Hook::PostIntegrationHook()
	{
		snapshot_->Changed(false);

		for(auto& cv : cvmanager_->GetCVs({}))
			cv->Evaluate(*snapshot_);

		for(auto& listener : listeners_)
			if(snapshot_->GetIteration() % listener->GetFrequency() == 0)
				listener->PostIntegration(snapshot_, cvmanager_);

		if(snapshot_->HasChanged())
			SyncToEngine();

		snapshot_->Changed(false);		
	}

	void Hook::PostSimulationHook()
	{
		snapshot_->Changed(false);

		for(auto& cv : cvmanager_->GetCVs({}))
			cv->Evaluate(*snapshot_);
		
		for(auto& listener : listeners_)
			listener->PostSimulation(snapshot_, cvs_);

		if(snapshot_->HasChanged())
			SyncToEngine();

		snapshot_->Changed(false);		
	}

	//! Sets the active snapshot.
	void Hook::SetSnapshot(Snapshot* snapshot)
	{
		snapshot_ = snapshot;
	}

	//! Sets the active Driver
	void Hook::SetDriver(Driver* driver)
	{
		driver_ = driver;
	}

	//! Sets the active CV manager.
	void Hook::SetCVManager(CVManager* cvmanager)
	{
		cvmanager_ = cvmanager;
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

	void Hook::NotifyObservers()
	{
		auto& comm = snapshot_->GetCommunicator();
		if(comm.rank() == 0)
			driver_->NotifyObservers(SimEvent(driver_, snapshot_->GetIteration()));
	}
}
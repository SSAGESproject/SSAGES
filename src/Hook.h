/**
 * This file is part of
 * SSAGES - Suite for Advanced Generalized Ensemble Simulations
 *
 * Copyright 2016 Hythem Sidky <hsidky@nd.edu>
 *                Ben Sikora <bsikora906@gmail.com>
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

#include "Simulations/SimObservable.h"
#include "EventListener.h"
#include "CVs/CollectiveVariable.h"
#include "Snapshot.h"
#include <algorithm>

namespace SSAGES
{
	//! Base class for hooks into the simultion engines
	/*!
	 * \ingroup core
	 *
	 * Abstract base class responsible for hooking into simulation engine
	 * and calling appropriate events.
	 */
	class Hook
	{
	private:
		//! Vector of event listeners.
		std::vector<EventListener*> _listeners;

		//! Vector of CVs.
		CVList _cvs;

		//! Driver running this hook
		Driver* _MDDriver;

	protected:
		//! Local snapshot.
		Snapshot* _snapshot;

		//! Synchronization to the simulation engine
		/*!
		 * A Hook must implement this method. It takes data from the snapshot
		 * and updates the simulation engine with it.
		 */
		virtual void SyncToEngine() = 0;

		//! Synchronization to the snapshot
		/*!
		 * A Hook must implement this method. It takes data from the simulation
		 * eingine and updates the snapshot with it.
		 */
		virtual void SyncToSnapshot() = 0;

	public:
		//! Constructor
		/*!
		 * Initialize a hook with world and walker communicators and
		 * corresponding walker ID.
		 */
		Hook() : 
		_listeners(0), _MDDriver(nullptr), _snapshot(nullptr)
		{}

		//! Sets the active snapshot.
		void SetSnapshot(Snapshot* snapshot);

		//! Sets the active Driver
		void SetMDDriver(Driver* MDDriver);

		//! Add a listener to the hook.
		/*!
		 * \param listener Pointer to the EventListener to be added to the Hook.
		 *
		 * Does nothing if the listener is already added.
		 */
		void AddListener(EventListener* listener);

		//! Add a CollectiveVariable to the hook.
		/*!
		 * \param cv CollectiveVariable to be added to the hook.
		 *
		 * Does nothing if the CollectiveVariable is already added.
		 */
		void AddCV(CollectiveVariable* cv);

		//! Notify observers of changes in the simulation.
		void NotifyObservers();
		//! Pre-simulation hook
		/*!
		 * This should be called at the appropriate
		 * time by the Hook implementation.
		 */
		void PreSimulationHook()
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

		//! Post-integration hook.
		/*!
		 * This function should be called by the Hook implementation within the
		 * integration routine such that the forces, position, velocities,
		 * etc.. will be updated.
		 */
		void PostIntegrationHook()
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

		//! Post-step hook.
		/*!
		 * This function should be called by the Hook implementation after the
		 * integration routine such that the forces, position, velocities,
		 * etc.. are done being updated.
		 */
		void PostStepHook()
		{
			NotifyObservers();
		}

		//! Post-simulation hook.
		/*!
		 * This method should be called by the Hook implementation at the
		 * end of the simulation.
		 */
		void PostSimulationHook()
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

		//! Destructor
		virtual ~Hook(){}
	};
}
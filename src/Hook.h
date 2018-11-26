/**
 * This file is part of
 * SSAGES - Software Suite for Advanced General Ensemble Simulations
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

#include "EventListener.h"
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
		std::vector<EventListener*> listeners_;

		//! Collective variable manager.
		class CVManager* cvmanager_;

	protected:
		//! Local snapshot.
		class Snapshot* snapshot_;

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
		listeners_(0), snapshot_(nullptr)
		{}

		//! Sets the active snapshot.
		void SetSnapshot(class Snapshot* snapshot);

		//! Sets the current CV manager.
		/*!
		 * \param cvmanager Collective variable manager.
		 */
		void SetCVManager(class CVManager* cvmanager);

		//! Add a listener to the hook.
		/*!
		 * \param listener Pointer to the EventListener to be added to the Hook.
		 *
		 * Does nothing if the listener is already added.
		 */
		void AddListener(EventListener* listener);

		//! Notify observers of changes in the simulation.
		void NotifyObservers();
		
		//! Pre-simulation hook
		/*!
		 * This should be called at the appropriate
		 * time by the Hook implementation.
		 */
		void PreSimulationHook();

		//! Post-integration hook.
		/*!
		 * This function should be called by the Hook implementation within the
		 * integration routine such that the forces, position, velocities,
		 * etc.. will be updated.
		 */
		void PostIntegrationHook();

		//! Post-step hook.
		/*!
		 * This function should be called by the Hook implementation after the
		 * integration routine such that the forces, position, velocities,
		 * etc.. are done being updated.
		 */
		void PostStepHook()
		{
		}

		//! Post-simulation hook.
		/*!
		 * This method should be called by the Hook implementation at the
		 * end of the simulation.
		 */
		void PostSimulationHook();

		//! Destructor
		virtual ~Hook(){}
	};
}
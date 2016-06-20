#pragma once

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

	public:
		//! Constructor
		/*!
		 * Initialize a hook with world and walker communicators and
		 * corresponding walker ID.
		 */
		Hook() : 
		_listeners(0), _snapshot(nullptr)
		{}

		//! Sets the active snapshot.
		void SetSnapshot(Snapshot* snapshot)
		{
			_snapshot = snapshot;
		}

		//! Add a listener to the hook.
		/*!
		 * \param listener Pointer to the EventListener to be added to the Hook.
		 *
		 * Does nothing if the listener is already added.
		 */
		void AddListener(EventListener* listener)
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
		void AddCV(CollectiveVariable* cv)
		{
			if(std::find(_cvs.begin(), _cvs.end(), cv) == _cvs.end())
				_cvs.push_back(cv);
		}

		//! Destructor
		virtual ~Hook(){}
	};
}
#pragma once

#include "EventListener.h"
#include "CVs/CollectiveVariable.h"
#include "Snapshot.h"
#include <algorithm>

namespace SSAGES
{
	// Abstract base class responsible for hooking into simulation engine 
	// and calling appropriate events.
	class Hook
	{
	private:
		// Vector of event listeners.
		std::vector<EventListener*> _listeners;

		// Vector of CVs.
		CVList _cvs;

	protected:
		// Local snapshot.
		Snapshot _snapshot;

		// Syncronization to engine. A Hook must implement this method 
		// where data is taken from the snapshot and updated within 
		// the simulation engine.
		virtual void SyncToEngine() = 0;

		// Synchronization to snapshot. A hook must implement this method 
		// where data is taken from the simulation eingine and updated 
		// within the snapshot.
		virtual void SyncToSnapshot() = 0;

		// Pre-simulation hook. This should be called at the appropriate 
		// time by the Hook implementation.
		void PreSimulationHook()
		{
			_snapshot.Changed(false);
		
			// Initialize/evaluate CVs.
			for(auto& cv : _cvs)
			{
				cv->Initialize(_snapshot);
				cv->Evaluate(_snapshot);
			}

			// Call presimulation method on listeners. 
			for(auto& listener : _listeners)
				listener->PreSimulation(&_snapshot, _cvs);

			// Sync snapshot to engine.
			if(_snapshot.HasChanged())
				SyncToEngine();

			_snapshot.Changed(false);
		}

		// Post-integration hook. This hould be called by the Hook 
		// implementation within the integration routine such that 
		// the forces, position, velocities, etc.. can be updated.
		void PostIntegrationHook()
		{
			_snapshot.Changed(false);

			for(auto& cv : _cvs)
				cv->Evaluate(_snapshot);

			for(auto& listener : _listeners){
				if(_snapshot.GetIteration() % listener->GetFrequency() == 0)
					listener->PostIntegration(&_snapshot, _cvs);
			}

			if(_snapshot.HasChanged())
				SyncToEngine();

			_snapshot.Changed(false);
		}

		// Post-simulation hook. This should be called by the Hook
		// implementation at the appropriate time.
		void PostSimulationHook()
		{
			_snapshot.Changed(false);

			for(auto& cv : _cvs)
				cv->Evaluate(_snapshot);
			
			for(auto& listener : _listeners)
				listener->PostSimulation(&_snapshot, _cvs);

			if(_snapshot.HasChanged())
				SyncToEngine();

			_snapshot.Changed(false);
		}

	public:
		// Initialize a hook.
		Hook() : 
		_listeners(0), _snapshot() {}

		// Add a listener to the hook. Does nothing 
		// if the listener is already added.
		void AddListener(EventListener* listener)
		{
			if(std::find(_listeners.begin(), _listeners.end(), listener) == _listeners.end())
				_listeners.push_back(listener);
		}

		// Adds a CV to the hook. Does nothing if the CV 
		// is already added.
		void AddCV(CollectiveVariable* cv)
		{
			if(std::find(_cvs.begin(), _cvs.end(), cv) == _cvs.end())
				_cvs.push_back(cv);
		}

		~Hook(){}
	};
}
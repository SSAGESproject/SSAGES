#pragma once

#include "EventListener.h"
#include "Snapshot.h"
#include <algorithm>

namespace SSAGES
{
	// Base class responsible for hooking into simulation engine 
	// and calling methods duing events.
	class Hook
	{
	private:
		std::vector<EventListener*> _listeners;

	protected:
		Snapshot _snapshot;

		virtual void SyncToEngine() = 0;

		virtual void SyncToSnapshot() = 0;

	public:
		Hook() : 
		_listeners(0), _snapshot() {}

		void PreSimulationHook()
		{
			for(auto& listener : _listeners)
				listener->PreSimulation(&_snapshot);

			if(_snapshot.HasChanged())
				SyncToEngine();

			_snapshot.Changed(false);
		}

		void PostIntegration()
		{
			for(auto& listener : _listeners)
				listener->PostIntegration(&_snapshot);

			if(_snapshot.HasChanged())
				SyncToEngine();

			_snapshot.Changed(false);
		}

		void PostSimulationHook()
		{
			for(auto& listener : _listeners)
				listener->PostSimulation(&_snapshot);

		}

		void AddListener(EventListener* listener)
		{
			if(std::find(_listeners.begin(), _listeners.end(), listener) == _listeners.end())
				_listeners.push_back(listener);
		}

		~Hook(){}
	};
}
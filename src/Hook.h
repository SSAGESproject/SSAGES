#pragma once

#include "EventListener.h"
#include "CollectiveVariable.h"
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
		CVList _cvs;

	protected:
		Snapshot _snapshot;

		virtual void SyncToEngine() = 0;

		virtual void SyncToSnapshot() = 0;

		void PreSimulationHook()
		{
			_snapshot.Changed(false);

			for(auto& listener : _listeners)
				listener->PreSimulation(&_snapshot, _cvs);

			if(_snapshot.HasChanged())
				SyncToEngine();

			_snapshot.Changed(false);
		}

		void PostIntegration()
		{
			_snapshot.Changed(false);

			for(auto& listener : _listeners){
				if(_snapshot.GetIteration() % listener->GetFrequency() == 0)
					listener->PostIntegration(&_snapshot, _cvs);
			}

			if(_snapshot.HasChanged())
				SyncToEngine();

			_snapshot.Changed(false);
		}

		void PostSimulationHook()
		{
			_snapshot.Changed(false);
			
			for(auto& listener : _listeners)
				listener->PostSimulation(&_snapshot, _cvs);

			if(_snapshot.HasChanged())
				SyncToEngine();

			_snapshot.Changed(false);
		}

	public:
		Hook() : 
		_listeners(0), _snapshot() {}

		void AddListener(EventListener* listener)
		{
			if(std::find(_listeners.begin(), _listeners.end(), listener) == _listeners.end())
				_listeners.push_back(listener);
		}

		void AddCV(CollectiveVariable* cv)
		{
			if(std::find(_cvs.begin(), _cvs.end(), cv) == _cvs.end())
				_cvs.push_back(cv);
		}

		~Hook(){}
	};
}
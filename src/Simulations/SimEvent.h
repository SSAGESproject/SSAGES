#pragma once

#include "../JSON/Serializable.h"

namespace SSAGES
{
	class SimObservable;

	class SimEvent
	{
		private:
			SimObservable* _observable;
			int _iteration;
			bool _forceObserve = false;

		public:
			SimEvent(SimObservable* observable, int iteration, bool forceObserve = false)
				: _observable(observable), _iteration(iteration), _forceObserve(forceObserve){}

			// Get Sim Observable.
			SimObservable const * GetObservable() const { return _observable; }

			// Get iteration.
			int GetIteration() const { return _iteration; }

			// Tell the observer that they should observe this event.
			bool ForceObserve() const { return _forceObserve; }
	};
}

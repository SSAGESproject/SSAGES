#pragma once

#include "../JSON/Serializable.h"

namespace SSAGES
{
	class SimObservable;

	//! Simulation event
	/*
	 * Simulation event that is passed to an observer.
	 *
	 * \ingroup Core
	 */
	class SimEvent
	{
		private:
			//! Pointer to the observable object from which the Event is originating.
			SimObservable* _observable;

			//! Iteration of the simulation.
			int _iteration;

			//! Force the observation of this event.
			bool _forceObserve = false;

		public:
			//! Constructor
			/*!
			 * \param observable Pointer to the observable object.
			 * \param iteration Interation.
			 * \param forceObserve If \c True force the observers to observe this event.
			 */
			SimEvent(SimObservable* observable, int iteration, bool forceObserve = false)
				: _observable(observable), _iteration(iteration), _forceObserve(forceObserve){}

			//! Get Sim Observable.
			/*!
			 * \return Pointer to the observable object.
			 */
			SimObservable const * GetObservable() const { return _observable; }

			//! Get iteration.
			/*!
			 * \return Iteration.
			 */
			int GetIteration() const { return _iteration; }

			//! Tell the observer that they should observe this event.
			/*!
			 * \return \c True if Observers are forced to observe the event.
			 */
			bool ForceObserve() const { return _forceObserve; }

			//! Destructor.
			virtual ~SimEvent();
	};
}

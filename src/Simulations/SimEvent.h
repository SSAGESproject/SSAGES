/**
 * This file has been adapted from
 * SAPHRON - Statistical Applied PHysics through Random On-the-fly Numerics
 * https://github.com/hsidky/SAPHRON
 *
 * Copyright 2016 Hythem Sidky
 *
 * This file is part of
 * SSAGES - Suite for Advanced Generalized Ensemble Simulations
 *
 * Copyright 2016 Ben Sikora <bsikora906@gmail.com>
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
			virtual ~SimEvent() {}
	};
}

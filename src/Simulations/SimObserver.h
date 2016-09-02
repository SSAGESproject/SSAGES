/**
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

#include "../Observers/Visitor.h"
#include "../JSON/Serializable.h"
#include "SimEvent.h"
#include <vector>
#include "json/json.h"
#include <boost/mpi.hpp>

namespace SSAGES
{
	class SimObserver; 

	//! List of Simulation Observers.
	typedef std::vector<SimObserver*> ObserverList;

	//! Base class for objects observing a simulation.
	/*
	 * Abstract class for objects wanting to observe a simulation.
	 *
	 * \ingroup Core
	 */
	class SimObserver : public Visitor, public Serializable
	{
		private:
			unsigned int _frequency = 1; //!< Frequency of observation.
			SimEvent _event; //!< Simulation event to be observed.

		protected:

			boost::mpi::communicator _world; //!< MPI global communicator.
			boost::mpi::communicator _comm; //!< MPI local communicator.
			int _numwalkers; //!< Number of walkers.
			int _wid; //!< ID of walker to be observed.

			//! Get current simulation event
			/*!
			 * \return Current simulation event.
			 */
			SimEvent& GetEvent()
			{
				return _event;
			}

			//! Set logging frequency.
			/*!
			 * \param f New frequency.
			 * \return New frequency.
			 */
			unsigned int SetFrequency(int f)
			{
				return _frequency = f;
			}

			//! Get current event iteration
			/*!
			 * \return Iteration of current event.
			 */
			unsigned int GetIteration()
			{
				return _event.GetIteration();
			}

			//! Get caller identifier.
			/*!
			 * \return 0 always.
			 */
			int GetObservableID()
			{
				return 0;
			}

		public:
			//! Constructor
			/*!
			 * \param world MPI global communicator.
			 * \param comm MPI local communicator.
			 * \param nws Number of walkers.
			 * \param wid ID of walker to observe.
			 * \param frequency Logging frequency.
			 *
			 * Initialize a SimObserver class with a specified observation
			 * frequency.
			 */
			SimObserver(boost::mpi::communicator world,
						boost::mpi::communicator comm,
						int nws,
						int wid,
						unsigned int frequency = 1)
				: _frequency(frequency), _event(nullptr, 0),
				_world(world), _comm(comm), _numwalkers(nws), _wid(wid){}

			//! Update observer when simulation has changed.
			/*!
			 * \param e New simulation event.
			 */
			void Update(const SimEvent& e);

			//! Get logging frequency.
			/*!
			 * \return Logging frequency.
			 */
			unsigned int GetFrequency() const { return _frequency; }

			//! Get name.
			/*!
			 * \return Name of the observer.
			 */
			virtual std::string GetName() const = 0;

			//! Called before visitors invokes.
			virtual void PreVisit(){}

			//! Called after visitors complete.
			virtual void PostVisit(){}

			//! Serialize observer.
			/*!
			 * \param json JSON value to store information to.
			 */
			virtual void Serialize(Json::Value& json) const = 0;

			//! Build the simulation observer.
			/*!
			 * \param json JSON value containing input information.
			 * \param world MPI global communicator.
			 * \param comm MPI local communicator.
			 * \param numwalkers Number of walkers.
			 * \param wid ID of the walker to observe.
			 * \param path Path for JSON path specification.
			 *
			 * \return Pointer to the newly created observer. \c nullptr in case
			 *         of unknown error.
			 *
			 * Static builder method for simobserver. If an unknown error
			 * occured, the return value is nullptr, but generally, this
			 * function will throw a BuildException on failure.
			 *
			 * \note Object lifetime is caller's responsibility!
			 */
			static SimObserver* BuildObserver(const Json::Value& json,
										boost::mpi::communicator world,
										boost::mpi::communicator comm,
										int numwalkers,
										int wid,
										const std::string& path);
			
			//! Build the simulation observer.
			/*!
			 * \param json JSON value containing input information.
			 * \param world MPI global communicator.
			 * \param comm MPI local communicator.
			 * \param numwalkers Number of walkers.
			 * \param wid ID of the walker to observe.
			 *
			 * \return Pointer to the newly created observer. \c nullptr in case
			 *         of unknown error.
			 *
			 * Static builder method for simobserver. If an unknown error
			 * occured, the return value is nullptr, but generally, this
			 * function will throw a BuildException on failure.
			 *
			 * \note Object lifetime is caller's responsibility!
			 */
			static SimObserver* BuildObserver(const Json::Value& json,
										boost::mpi::communicator world,
										boost::mpi::communicator comm,
										int numwalkers,
										int wid);
			
			//! Build list of observers.
			/*!
			 * \param json JSON value containing input information.
			 * \param world MPI global communicator.
			 * \param comm MPI local communicator.
			 * \param numwalkers Number of walkers.
			 * \param wid ID of the walker to be observed.
			 * \param ol List of observers to which the new observers are appended.
			 *
			 * Static builder method for sim observers. ObserverList is a vector
			 * which will be populated with initialized observers.
			 *
			 * \note Object lifetime is caller's responsibility!
			 */
			static void BuildObservers(const Json::Value& json,
			boost::mpi::communicator world,
			boost::mpi::communicator comm,
			int numwalkers,
			int wid,
			ObserverList& ol);

			//! Destructor.
			virtual ~SimObserver(){}	
	};
}

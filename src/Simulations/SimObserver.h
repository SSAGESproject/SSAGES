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
	typedef std::vector<SimObserver*> ObserverList;

	// Abstract class for objects wanting to observe a simulation.
	class SimObserver : public Visitor, public Serializable
	{
		private:
			unsigned int _frequency = 1;
			SimEvent _event;

		protected:

			boost::mpi::communicator _world;
			boost::mpi::communicator _comm;
			int _numwalkers;
			int _wid;

			SimEvent& GetEvent()
			{
				return _event;
			}

			// Set logging frequency.
			unsigned int SetFrequency(int f)
			{
				return _frequency = f;
			}

			// Get current event iteration
			unsigned int GetIteration()
			{
				return _event.GetIteration();
			}

			// Get caller identifier.
			int GetObservableID()
			{
				return 0;
			}

		public:
			// Initialize a SimObserver class with a specified observation frequency.
			SimObserver(boost::mpi::communicator world,
						boost::mpi::communicator comm,
						int nws,
						int wid,
						unsigned int frequency = 1)
				: _frequency(frequency), _event(nullptr, 0),
				_world(world), _comm(comm), _numwalkers(nws), _wid(wid){}

			// Update observer when simulation has changed.
			void Update(const SimEvent& e);


			// Get logging frequency.
			unsigned int GetFrequency() const { return _frequency; }

			// Get name. 
			virtual std::string GetName() const = 0;

			// Called before visitors invokes.
			virtual void PreVisit(){}

			// Called after visitors complete.
			virtual void PostVisit(){}

			// Serialize observer.
			virtual void Serialize(Json::Value& json) const = 0;

			// Static builder method for simobserver. If return value is nullptr, 
			// then an unknown error occurred. It will throw a BuildException on 
			// failure.  Object lifetime is caller's responsibility!
			static SimObserver* BuildObserver(const Json::Value& json,
										boost::mpi::communicator world,
										boost::mpi::communicator comm,
										int numwalkers,
										int wid,
										const std::string& path);
			
			// Static builder method for simobserver. If return value is nullptr, 
			// then an unknown error occurred. It will throw a BuildException on 
			// failure.  Object lifetime is caller's responsibility!
			static SimObserver* BuildObserver(const Json::Value& json,
										boost::mpi::communicator world,
										boost::mpi::communicator comm,
										int numwalkers,
										int wid);
			
			// Static builder method for sim observers. ObserverList is a vector
			// which will be populated with initialized observers. Object lifetime
			// is caller's responsibility!
			static void BuildObservers(const Json::Value& json,
			boost::mpi::communicator world,
			boost::mpi::communicator comm,
			int numwalkers,
			int wid,
			ObserverList& ol);

			virtual ~SimObserver(){}	
	};
}

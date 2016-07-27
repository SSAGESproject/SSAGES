#include "Hook.h"
#include "Drivers/Driver.h"
#include <algorithm>

namespace SSAGES
{
	//! Sets the active snapshot.
	void Hook::SetSnapshot(Snapshot* snapshot)
	{
		_snapshot = snapshot;
	}

	//! Sets the active Driver
	void Hook::SetMDDriver(Driver* MDDriver)
	{
		_MDDriver = MDDriver;
	}

	//! Add a listener to the hook.
	/*!
	 * \param listener Pointer to the EventListener to be added to the Hook.
	 *
	 * Does nothing if the listener is already added.
	 */
	void Hook::AddListener(EventListener* listener)
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
	void Hook::AddCV(CollectiveVariable* cv)
	{
		if(std::find(_cvs.begin(), _cvs.end(), cv) == _cvs.end())
			_cvs.push_back(cv);
	}

	void Hook::NotifyObservers()
	{
		auto& comm = _snapshot->GetCommunicator();
		if(comm.rank() == 0)
			_MDDriver->NotifyObservers(SimEvent(_MDDriver, _MDDriver->GetIteration()));
	}
}
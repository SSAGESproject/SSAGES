#pragma once 

#include "types.h"

namespace SSAGES
{
	// Forward declare. 
	class Snapshot;

	//! Base abstract class for listening in to events fired by "Hook".
	/*!
	 * \ingroup Core
	 */
	class EventListener
	{
	private:
		unsigned int _frequency; //!< Frequency for listening.
	
	public:
		//! Constructor
		/*!
		 * \param frequency Frequency for listening.
		 */
		EventListener(unsigned int frequency) : 
		_frequency(frequency)
		{

		}

		//! Get frequency of event listener.
		/*!
		 * \return Frequency of event listener.
		 */
		unsigned int GetFrequency() const { return _frequency; }

		//! Method call prior to simulation initiation.
		/*!
		 * \param snapshot Pointer to the simulation snapshot.
		 * \param cvs List of CVs.
		 *
		 * This function will be called before the simulation is started.
		 */
		virtual void PreSimulation(Snapshot* snapshot, const CVList& cvs) = 0;

		//! Method call post integration.
		/*!
		 * \param snapshot Pointer to the simulation snapshot.
		 * \param cvs List of CVs.
		 *
		 * This function will be called at the end of each integration step.
		 */
		virtual void PostIntegration(Snapshot* snapshot, const CVList& cvs) = 0;

		//! Method call post simulation.
		/*!
		 * \param snapshot Pointer to the simulation snapshot.
		 * \param cvs List of CVs.
		 *
		 * This function will be called after the simulation has finished.
		 */
		virtual void PostSimulation(Snapshot* snapshot, const CVList& cvs) = 0;

		//! Destructor
		virtual ~EventListener() { }
	};
}
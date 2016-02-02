#ifdef FIX_CLASS

FixStyle(ssages,FixSSAGES)

#else
#ifndef LMP_FIX_SSAGES_H
#define LMP_FIX_SSAGES_H

#include "fix.h"
#include "Hook.h"

namespace LAMMPS_NS
{
	// SSAGES Hook class for LAMMPS implemented as 
	// a LAMMPS fix. This is activated by adding 
	// a "ssages" fix to "all". Note that thermo must 
	// be set to 1 in order for the synchronizing to work.
	class FixSSAGES : public Fix, SSAGES::Hook
	{
	protected:
		// Implementation of the SyncToEngine interface.
		void SyncToEngine() override;

		// Implementation of the SyncToSnapshot interface.
		void SyncToSnapshot() override;

	public:
		FixSSAGES(class LAMMPS *, int, char**);

		// Setup for presimulation call.
		void setup(int) override;
  		
		// Post force where the synchronization occurs.
  		void post_force(int) override;

  		// Post-run for post-simulation call.
  		void post_run() override;

  		// Set mask to let LAMMPS know what triggers we're interested in.
  		int setmask() override;
	};
}

#endif
#endif
#ifdef FIX_CLASS

FixStyle(ssages,FixSSAGES)

#else
#ifndef LMP_FIX_SSAGES_H
#define LMP_FIX_SSAGES_H

#include "fix.h"
#include "Hook.h"


namespace LAMMPS_NS
{
  inline const std::array<double, 6> ConvertToLatticeConstant(const std::array<double, 6>& lmmpsbox)
  {
   
    std::array<double, 6> box;

    box[0] = lmmpsbox[0];   
    box[1] = sqrt(lmmpsbox[1]*lmmpsbox[1] + lmmpsbox[3]*lmmpsbox[3]);
    box[2] = sqrt(lmmpsbox[2]*lmmpsbox[2] + lmmpsbox[4]*lmmpsbox[4] + lmmpsbox[5]*lmmpsbox[5]);
    box[3] = acos((lmmpsbox[3]*lmmpsbox[4] + lmmpsbox[1]*lmmpsbox[5])/box[1]*box[2]);
    box[4] = acos(lmmpsbox[4]/box[2]);
    box[5] = acos(lmmpsbox[3]/box[1]);

    return box;
  }
	// SSAGES Hook class for LAMMPS implemented as 
	// a LAMMPS fix. This is activated by adding 
	// a "ssages" fix to "all". Note that thermo must 
	// be set to 1 in order for the synchronizing to work.
	class FixSSAGES : public Fix, public SSAGES::Hook
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
  		
  		// Post-step for post-step call.
  		void end_of_step() override;

  		// Set mask to let LAMMPS know what triggers we're interested in.
  		int setmask() override;

  		// Grab the box info from lammps
  		const std::array<double, 6> GatherLAMMPSVectors() const;
	};
}

#endif
#endif
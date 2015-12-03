#include "fix_ssages.h"
#include "atom.h"
#include <iostream>

using namespace SSAGES;
using namespace LAMMPS_NS::FixConst;

namespace LAMMPS_NS
{
	FixSSAGES::FixSSAGES(LAMMPS *lmp, int narg, char **arg) : 
	Fix(lmp, narg, arg), Hook()
	{
	}

	void FixSSAGES::setup(int)
	{
		// Allocate vectors for snapshot.
		auto n = atom->natoms;

		SyncToSnapshot();
		SyncToEngine();
		Hook::PreSimulationHook();
	}

	int FixSSAGES::setmask()
	{
	  int mask = 0;
	  mask |= PRE_EXCHANGE;
	  return mask;
	}

	void FixSSAGES::SyncToSnapshot()
	{

	}

	void FixSSAGES::SyncToEngine()
	{

	}

	void FixSSAGES::post_force(int)
	{

	}
}
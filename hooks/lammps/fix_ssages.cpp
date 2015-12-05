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
		auto n = atom->natoms + 1;
		auto& pos = _snapshot.GetPositions();
		pos.resize(n);
		auto& vel = _snapshot.GetVelocities();
		vel.resize(n);
		auto& frc = _snapshot.GetForces();
		frc.resize(n);
		auto& ids = _snapshot.GetAtomIDs();
		ids.resize(n);
		auto& mids = _snapshot.GetMoleculeIDs();
		mids.resize(n);
		auto& types = _snapshot.GetAtomTypes();
		types.resize(n);

		SyncToSnapshot();
		Hook::PreSimulationHook();
	}

	void FixSSAGES::post_force(int)
	{
		SyncToSnapshot();
		Hook::PostIntegration();
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
}
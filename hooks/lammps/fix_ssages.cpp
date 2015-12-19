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

	//Will force a call to SyncToEngine?
	void FixSSAGES::SyncToSnapshot() //put LAMMPS values -> Snapshot
	{

		// Obtain local reference to snapshot variables.
		// Const Atom will ensure that atom variables are
		// not being changed. Only snapshot side variables should
		// change.
		const Atom& _atom=atom;

		auto& pos = _snapshot.GetPositions();
		auto& vel = _snapshot.GetVelocities();
		auto& frc = _snapshot.GetForces();

		// Labels and ids for future work on only updating
		// atoms that have changed.
		auto& ids = _snapshot.GetAtomIDs();
		auto& mids = _snapshot.GetMoleculeIDs();
		auto& types = _snapshot.GetAtomTypes();

		// Positions
		for (int i=0; i<_atom->x.size(); i++)
		{
			for (int j=0; j<_atom->x[i].size(); j++)
			{
				pos[3*i+j]=_atom->x[i][j];
			}
		}

		// Forces
		for (int i=0; i<_atom->f.size(); i++)
		{
			for (int j=0; j<_atom->f[i].size(); j++)
			{
				frc[3*i+j]=_atom->f[i][j];
			}
		}

		//velocities
		for (int i=0; i<_atom->v.size(); i++)
		{
			for (int j=0; j<_atom->v[i].size(); j++)
			{
				vel[3*i+j]=_atom->v[i][j];
			}
		}
	}

	void FixSSAGES::SyncToEngine() //put Snapshot values -> LAMMPS
	{
		// Obtain local const reference to snapshot variables.
		// Const will ensure that _snapshot variables are
		// not being changed. Only engine side variables should
		// change. 

		const Vector& pos = _snapshot.GetPositions();
		const Vector& vel = _snapshot.GetVelocities();
		const Vector& frc = _snapshot.GetForces();

		// Labels and ids for future work on only updating
		// atoms that have changed.
		const Label& ids = _snapshot.GetAtomIDs();
		const Label& mids = _snapshot.GetMoleculeIDs();
		const Label& types = _snapshot.GetAtomTypes();

		// Need to parallize this portion each processor can set and read
		// certain atoms.
		// Loop through all atoms and set their values
		//Positions
		for (int i=0; i<atom->x.size(); i++)
		{
			for (int j=0; j<atom->x[i].size(); j++)
			{
				atom->x[i][j]=pos[3*i+j];
			}
		}

		// Forces
		for (int i=0; i<atom->f.size(); i++)
		{
			for (int j=0; j<atom->f[i].size(); j++)
			{
				atom->f[i][j]=frc[3*i+j];
			}
		}

		//velocities
		for (int i=0; i<atom->v.size(); i++)
		{
			for (int j=0; j<atom->v[i].size(); j++)
			{
				atom->v[i][j]=vel[3*i+j];
			}
		}

		// Set thermodynamic information

	}
}
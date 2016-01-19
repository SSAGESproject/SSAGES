#include <iostream>
#include "fix_ssages.h"
#include "atom.h"
#include "compute.h"
#include "modify.h"
#include "force.h"
#include "update.h"
#include "domain.h"
#include "MockMethod.h"

using namespace SSAGES;
using namespace LAMMPS_NS::FixConst;

namespace LAMMPS_NS
{
	FixSSAGES::FixSSAGES(LAMMPS *lmp, int narg, char **arg) : 
	Fix(lmp, narg, arg), Hook()
	{
		this->AddListener(new MockMethod(1));
	}

	void FixSSAGES::setup(int)
	{
		// Allocate vectors for snapshot.
		auto n = atom->natoms;
		auto& pos = _snapshot.GetPositions();
		pos.resize(n);
		auto& vel = _snapshot.GetVelocities();
		vel.resize(n);
		auto& frc = _snapshot.GetForces();
		frc.resize(n);
		auto& ids = _snapshot.GetAtomIDs();
		ids.resize(n);
		auto& types = _snapshot.GetAtomTypes();
		types.resize(n);

		SyncToSnapshot();
		Hook::PreSimulationHook();
	}

	void FixSSAGES::pre_force(int)
	{
		SyncToSnapshot();
		Hook::PostIntegration();
	}

	int FixSSAGES::setmask()
	{
	  int mask = 0;
	  mask |= PRE_FORCE;
	  return mask;
	}

	void FixSSAGES::SyncToSnapshot() //put LAMMPS values -> Snapshot
	{

		// Obtain local reference to snapshot variables.
		// Const Atom will ensure that atom variables are
		// not being changed. Only snapshot side variables should
		// change.
		const auto* _atom = atom;

		auto& pos = _snapshot.GetPositions();
		auto& vel = _snapshot.GetVelocities();
		auto& frc = _snapshot.GetForces();

		// Labels and ids for future work on only updating
		// atoms that have changed.
		auto& ids = _snapshot.GetAtomIDs();
		auto& types = _snapshot.GetAtomTypes();

		// Thermo properties:
		Compute *thermoproperty;
		int icompute; 

		//Temperature
		const char* id_temp = "thermo_temp";
		icompute = modify->find_compute(id_temp);
		auto* temperature = modify->compute[icompute];
		_snapshot.SetTemperature(temperature->compute_scalar());
		
		//Pressure
		const char* id_press = "thermo_press";
		icompute = modify->find_compute(id_press);
		thermoproperty = modify->compute[icompute];
		_snapshot.SetPressure(thermoproperty->compute_scalar());
		
		//Energy
		double etot = 0;

		// Get potential energy.
		const char* id_pe = "thermo_pe";
		icompute = modify->find_compute(id_pe);
		auto* pe = modify->compute[icompute];
		etot += pe->scalar;

		// Compute kinetic energy.
		double ekin = 0.5 * temperature->scalar * temperature->dof  * force->boltz;
		etot += ekin;

		// Store in snapshot.
		_snapshot.SetEnergy(etot/atom->natoms);
		
		// Get iteration.
		_snapshot.SetIteration(update->ntimestep);
		
		// Get volume.
		double vol = 0.;
		if (domain->dimension == 3)
			vol = domain->xprd * domain->yprd * domain->zprd;
		else
			vol = domain->xprd * domain->yprd;
		_snapshot.SetVolume(vol);

		// Positions
		for (int i = 0; i < 3; ++i)
		{
			pos[i][0] = _atom->x[i][0]; //x
			pos[i][1] = _atom->x[i][1]; //y
			pos[i][2] = _atom->x[i][2]; //z
			
			frc[i][0] = _atom->f[i][0]; //force->x
			frc[i][1] = _atom->f[i][1]; //force->y
			frc[i][2] = _atom->f[i][2]; //force->z
			
			vel[i][0] = _atom->v[i][0];
			vel[i][1] = _atom->v[i][1];
			vel[i][2] = _atom->v[i][2];
			
			ids[i] = _atom->tag[i];
			
			types[i] = _atom->type[i];
		}
	}

	void FixSSAGES::SyncToEngine() //put Snapshot values -> LAMMPS
	{
		// Obtain local const reference to snapshot variables.
		// Const will ensure that _snapshot variables are
		// not being changed. Only engine side variables should
		// change. 
		const auto& pos = _snapshot.GetPositions();
		const auto& vel = _snapshot.GetVelocities();
		const auto& frc = _snapshot.GetForces();

		// Labels and ids for future work on only updating
		// atoms that have changed.
		const auto& ids = _snapshot.GetAtomIDs();
		const auto& types = _snapshot.GetAtomTypes();

		// Loop through all atoms and set their values
		// Positions
		for (int i = 0; i < atom->natoms; ++i)
		{
			atom->x[i][0] = pos[i][0]; //x 
			atom->x[i][1] = pos[i][1]; //y
			atom->x[i][2] = pos[i][2]; //z
			atom->f[i][0] = frc[i][0]; //force->x
			atom->f[i][1] = frc[i][1]; //force->y
			atom->f[i][2] = frc[i][2]; //force->z
			atom->v[i][0] = vel[i][0]; //velocity->x
			atom->v[i][1] = vel[i][1]; //velocity->y
			atom->v[i][2] = vel[i][2]; //velocity->z
			atom->tag[i] = ids[i];
			atom->type[i] = types[i];
		}

		// LAMMPS computes will reset thermo data based on
		// updated information. No need to sync thermo data
		// from snapshot to engine.
		// However, this could change in the future.
	}
}
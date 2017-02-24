#include <iostream>
#include <numeric>
#include "fix_ssages.h"
#include "atom.h"
#include "compute.h"
#include "modify.h"
#include "force.h"
#include "pair.h"
#include "update.h"
#include "domain.h"
#include <boost/mpi.hpp>
#include <boost/version.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>

using namespace SSAGES;
using namespace LAMMPS_NS::FixConst;
using namespace boost;

namespace LAMMPS_NS
{
	FixSSAGES::FixSSAGES(LAMMPS *lmp, int narg, char **arg) : 
	Fix(lmp, narg, arg), Hook()
	{
		const char* idtemp = "thermo_temp";
		int itemp = modify->find_compute(idtemp);
		tempid_ = modify->compute[itemp];
		tempid_->addstep(update->ntimestep + 1);

		const char* idpe = "thermo_pe";
		int ipe = modify->find_compute(idpe);
		peid_ = modify->compute[ipe];
		peid_->addstep(update->ntimestep + 1);

		// Turn on virial flag. 
		virial_flag = 1;
	}

	void FixSSAGES::setup(int)
	{
		// Allocate vectors for snapshot.
		// Here we size to local number of atoms. We will 
		// resize later.
		auto n = atom->nlocal;

		auto& pos = snapshot_->GetPositions();
		pos.resize(n);
		auto& vel = snapshot_->GetVelocities();
		vel.resize(n);
		auto& frc = snapshot_->GetForces();
		frc.resize(n);
		auto& masses = snapshot_->GetMasses();
		masses.resize(n);
		auto& ids = snapshot_->GetAtomIDs();
		ids.resize(n);
		auto& types = snapshot_->GetAtomTypes();
		types.resize(n);
		auto& charges = snapshot_->GetCharges();
		charges.resize(n);

		auto& flags = snapshot_->GetImageFlags();
		flags.resize(n);

		SyncToSnapshot();
		Hook::PreSimulationHook();
	}

	void FixSSAGES::post_force(int)
	{
		SyncToSnapshot();
		Hook::PostIntegrationHook();
	}

	void FixSSAGES::post_run()
	{
		SyncToSnapshot();
		Hook::PostSimulationHook();
	}

	void FixSSAGES::end_of_step()
	{
		Hook::PostStepHook();
	}

	int FixSSAGES::setmask()
	{
		// We are interested in post-force and post run hooks.
		int mask = 0;
		mask |= POST_FORCE;
		mask |= POST_RUN;
		mask |= END_OF_STEP;
		return mask;
	}

	void FixSSAGES::SyncToSnapshot() //put LAMMPS values -> Snapshot
	{
		using namespace boost::mpi;

		// Obtain local reference to snapshot variables.
		// Const Atom will ensure that atom variables are
		// not being changed. Only snapshot side variables should
		// change.
		const auto* atom_ = atom;

		auto n = atom->nlocal;

		auto& pos = snapshot_->GetPositions();
		pos.resize(n);
		auto& vel = snapshot_->GetVelocities();
		vel.resize(n);
		auto& frc = snapshot_->GetForces();
		frc.resize(n);
		auto& masses = snapshot_->GetMasses();
		masses.resize(n);
		auto& flags = snapshot_->GetImageFlags();
		flags.resize(n);

		// Labels and ids for future work on only updating
		// atoms that have changed.
		auto& ids = snapshot_->GetAtomIDs();
		ids.resize(n);
		auto& types = snapshot_->GetAtomTypes();
		types.resize(n);

		auto& charges = snapshot_->GetCharges();
		charges.resize(n);

		// Thermo properties:
		snapshot_->SetTemperature(tempid_->compute_scalar());
		tempid_->addstep(update->ntimestep + 1);
		
		//Energy
		double etot = 0;

		// Get potential energy.
		etot += peid_->compute_scalar();
		peid_->addstep(update->ntimestep + 1);

		// Compute kinetic energy.
		double ekin = 0.5 * tempid_->scalar * tempid_->dof  * force->boltz;
		etot += ekin;

		// Store in snapshot.
		snapshot_->SetEnergy(etot/atom->natoms);
		snapshot_->SetIteration(update->ntimestep);
		snapshot_->SetDielectric(force->dielectric);
		snapshot_->Setqqrd2e(force->qqrd2e);
		snapshot_->SetNumAtoms(n);
		
		// Get H-matrix.
		Matrix3 H;
		H << domain->h[0], domain->h[5], domain->h[4],
		                0, domain->h[1], domain->h[3],
		                0,            0, domain->h[2];

		snapshot_->SetHMatrix(H);
		snapshot_->SetKb(force->boltz);

		// Get box origin. 
		Vector3 origin;
		origin = {
			domain->boxlo[0], 
			domain->boxlo[1], 
			domain->boxlo[2]
		};
	
		snapshot_->SetOrigin(origin);

		// Set periodicity. 
		snapshot_->SetPeriodicity({
			domain->xperiodic, 
			domain->yperiodic, 
			domain->zperiodic
		});

		// Zero the virial - we are only interested in accumulation.
		snapshot_->SetVirial(Matrix3::Zero());

		// First we sync local data, then gather.
		// we gather data across all processors.
		for (int i = 0; i < n; ++i)
		{
			pos[i][0] = atom_->x[i][0]; //x
			pos[i][1] = atom_->x[i][1]; //y
			pos[i][2] = atom_->x[i][2]; //z
			
			frc[i][0] = atom_->f[i][0]; //force->x
			frc[i][1] = atom_->f[i][1]; //force->y
			frc[i][2] = atom_->f[i][2]; //force->z

			if(atom_->rmass_flag)
				masses[i] = atom_->rmass[i];
			else
				masses[i] = atom_->mass[atom_->type[i]];
			
			vel[i][0] = atom_->v[i][0];
			vel[i][1] = atom_->v[i][1];
			vel[i][2] = atom_->v[i][2];

			// Image flags. 
			flags[i][0] = (atom_->image[i] & IMGMASK) - IMGMAX;;
			flags[i][1] = (atom_->image[i] >> IMGBITS & IMGMASK) - IMGMAX;
			flags[i][2] = (atom_->image[i] >> IMG2BITS) - IMGMAX;

			ids[i] = atom_->tag[i];
			types[i] = atom_->type[i];

			if(atom_->q_flag)
				charges[i] = atom_->q[i];
			else
				charges[i] = 0;
		}
	}

	void FixSSAGES::SyncToEngine() //put Snapshot values -> LAMMPS
	{
		// Obtain local const reference to snapshot variables.
		// Const will ensure that snapshot_ variables are
		// not being changed. Only engine side variables should
		// change. 
		const auto& pos = snapshot_->GetPositions();
		const auto& vel = snapshot_->GetVelocities();
		const auto& frc = snapshot_->GetForces();
		const auto& masses = snapshot_->GetMasses();
		const auto& charges = snapshot_->GetCharges();

		// Labels and ids for future work on only updating
		// atoms that have changed.
		const auto& ids = snapshot_->GetAtomIDs();
		const auto& types = snapshot_->GetAtomTypes();

		for (int i = 0; i < atom->nlocal; ++i)
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
			
			if(atom->q_flag)
				atom->q[i] = charges[i];
			
			// Current masses can only be changed if using
			// per atom masses in lammps
			if(atom->rmass_flag)
				atom->rmass[i] = masses[i];			
		}

		// Update the virial (negative contribution to box). 
		const auto& vir = snapshot_->GetVirial();
		virial[0] = -vir(0,0);
		virial[1] = -vir(1,1);
		virial[2] = -vir(2,2);
		virial[3] = -vir(0,1);
		virial[4] = -vir(0,2);
		virial[5] = -vir(1,2);
	}

	FixSSAGES::~FixSSAGES()
	{
	}
}


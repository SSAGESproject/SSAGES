#include <iostream>
#include "fix_ssages.h"
#include "atom.h"
#include "compute.h"
#include "modify.h"
#include "force.h"
#include "update.h"
#include "domain.h"
#include <boost/mpi.hpp>

using namespace SSAGES;
using namespace LAMMPS_NS::FixConst;
using namespace boost;

namespace LAMMPS_NS
{
	// Copyright (C) 2015 Lorenz Hübschle-Schneider <lorenz@4z2.de>
	// Helper function for variable MPI_allgather.
	template <typename T>
	void allgatherv_serialize(const mpi::communicator &comm, const std::vector<T> &in, std::vector<T> &out) 
	{
		// Step 1: serialize input data
		mpi::packed_oarchive oa(comm);
		oa << in;

		// Step 2: exchange sizes (archives' .size() is measured in bytes)
		// Need to cast to int because this is what MPI uses as size_t...
		const int in_size = static_cast<int>(in.size()),
				transmit_size = static_cast<int>(oa.size());

		std::vector<int> in_sizes(comm.size()), transmit_sizes(comm.size());
		mpi::all_gather(comm, in_size, in_sizes.data());
		mpi::all_gather(comm, transmit_size, transmit_sizes.data());

		// Step 3: calculate displacements from sizes (prefix sum)
		std::vector<int> displacements(comm.size() + 1);
		displacements[0] = sizeof(mpi::packed_iarchive);
		for (int i = 1; i <= comm.size(); ++i)
			displacements[i] = displacements[i-1] + transmit_sizes[i-1];
		

		// Step 4: allocate space for result and MPI_Allgatherv
		char* recv = new char[displacements.back()];
		auto sendptr = const_cast<void*>(oa.address());
		auto sendsize = oa.size();

		int status = MPI_Allgatherv(sendptr, sendsize, MPI_PACKED, recv,
		                            transmit_sizes.data(), displacements.data(),
		                            MPI_PACKED, comm);

		if (status != 0) 
		{
		    std::cerr << "PE " << comm.rank() << ": MPI_Allgatherv returned "
		              << status << ", errno " << errno << std::endl;
		    return;
		}

		// Step 5: deserialize received data
		// Preallocate storage to prevent reallocations
		std::vector<T> temp;
		size_t largest_size = *std::max_element(in_sizes.begin(), in_sizes.end());
		temp.reserve(largest_size);
		out.reserve(std::accumulate(in_sizes.begin(), in_sizes.end(), 0));

		// Deserialize archives one by one, inserting elements at the end.
		for (int i = 0; i < comm.size(); ++i) 
		{
			mpi::packed_iarchive archive(comm);
			archive.resize(transmit_sizes[i]);
			memcpy(archive.address(), recv + displacements[i], transmit_sizes[i]);

			temp.clear();
			temp.resize(in_sizes[i]);
			archive >> temp;
			out.insert(out.end(), temp.begin(), temp.end());
		}
    }

	FixSSAGES::FixSSAGES(LAMMPS *lmp, int narg, char **arg) : 
	Fix(lmp, narg, arg), Hook()
	{
	}

	void FixSSAGES::setup(int)
	{
		// Allocate vectors for snapshot.
		// Here we size to local number of atoms. We will 
		// resize later.
		auto n = atom->nlocal;
		auto& pos = _snapshot->GetPositions();
		pos.resize(n);
		auto& vel = _snapshot->GetVelocities();
		vel.resize(n);
		auto& frc = _snapshot->GetForces();
		frc.resize(n);
		auto& ids = _snapshot->GetAtomIDs();
		ids.resize(n);
		auto& types = _snapshot->GetAtomTypes();
		types.resize(n);

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

	int FixSSAGES::setmask()
	{
		// We are interested in post-force and post run hooks.
		int mask = 0;
		mask |= POST_FORCE;
		mask |= POST_RUN;
		return mask;
	}

	void FixSSAGES::SyncToSnapshot() //put LAMMPS values -> Snapshot
	{
		using namespace boost::mpi;

		// Obtain local reference to snapshot variables.
		// Const Atom will ensure that atom variables are
		// not being changed. Only snapshot side variables should
		// change.
		const auto* _atom = atom;

		auto n = atom->nlocal;

		auto& pos = _snapshot->GetPositions();
		pos.resize(n);
		auto& vel = _snapshot->GetVelocities();
		vel.resize(n);
		auto& frc = _snapshot->GetForces();
		frc.resize(n);

		// Labels and ids for future work on only updating
		// atoms that have changed.
		auto& ids = _snapshot->GetAtomIDs();
		ids.resize(n);
		auto& types = _snapshot->GetAtomTypes();
		types.resize(n);

		// Thermo properties:
		Compute *thermoproperty;
		int icompute; 

		//Temperature
		const char* id_temp = "thermo_temp";
		icompute = modify->find_compute(id_temp);
		auto* temperature = modify->compute[icompute];
		_snapshot->SetTemperature(temperature->compute_scalar());
		
		//Pressure
		const char* id_press = "thermo_press";
		icompute = modify->find_compute(id_press);
		thermoproperty = modify->compute[icompute];
		_snapshot->SetPressure(thermoproperty->compute_scalar());
		
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
		_snapshot->SetEnergy(etot/atom->natoms);
		
		// Get iteration.
		_snapshot->SetIteration(update->ntimestep);
		
		// Get volume.
		double vol = 0.;
		if (domain->dimension == 3)
			vol = domain->xprd * domain->yprd * domain->zprd;
		else
			vol = domain->xprd * domain->yprd;
		_snapshot->SetVolume(vol);

		// First we sync local data, then gather.
		// we gather data across all processors.
		for (int i = 0; i < n; ++i)
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

		auto& comm = _snapshot->GetCommunicator();
		allgatherv_serialize(comm, pos, pos);
		allgatherv_serialize(comm, vel, vel);
		allgatherv_serialize(comm, frc, frc);
		allgatherv_serialize(comm, ids, ids);
		allgatherv_serialize(comm, types, types);
	}

	void FixSSAGES::SyncToEngine() //put Snapshot values -> LAMMPS
	{
		// Obtain local const reference to snapshot variables.
		// Const will ensure that _snapshot variables are
		// not being changed. Only engine side variables should
		// change. 
		const auto& pos = _snapshot->GetPositions();
		const auto& vel = _snapshot->GetVelocities();
		const auto& frc = _snapshot->GetForces();

		// Labels and ids for future work on only updating
		// atoms that have changed.
		const auto& ids = _snapshot->GetAtomIDs();
		const auto& types = _snapshot->GetAtomTypes();

		// Loop through all atoms and set their values
		// Positions. This is predicatd on the fact that the 
		// allgatherv_serialize keeps local processor information 
		// FIRST. This means that the first nlocal entries are the 
		// data corresponding to that on the local processor.
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
		}

		// LAMMPS computes will reset thermo data based on
		// updated information. No need to sync thermo data
		// from snapshot to engine.
		// However, this will change in the future.
	}
}
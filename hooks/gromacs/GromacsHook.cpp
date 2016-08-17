#include "GromacsHook.h"
#include <boost/mpi.hpp>

using namespace boost;

namespace SSAGES
{
	// Copyright (C) 2015 Lorenz Hübschle-Schneider <lorenz@4z2.de>
	// Helper function for variable MPI_allgather.
	template <typename T, typename transmit_type = uint64_t>
	void allgatherv_serialize(const mpi::communicator &comm, std::vector<T> &in, std::vector<T> &out) 
	{
		// Step 1: exchange sizes
		// We need to compute the displacement array, specifying for each PE
		// at which position in out to place the data received from it
		// Need to cast to int because this is what MPI uses as size_t...
		const int factor = sizeof(T) / sizeof(transmit_type);
		const int in_size = static_cast<int>(in.size()) * factor;
		std::vector<int> sizes(comm.size());
		mpi::all_gather(comm, in_size, sizes.data());

		// Step 2: calculate displacements from sizes
		// Compute prefix sum to compute displacements from sizes
		std::vector<int> displacements(comm.size() + 1);
		displacements[0] = 0;
		std::partial_sum(sizes.begin(), sizes.end(), displacements.begin() + 1);

		// divide by factor by which T is larger than transmit_type
		out.resize(displacements.back() / factor);
		
		// Step 3: MPI_Allgatherv
		transmit_type *sendptr = reinterpret_cast<transmit_type*>(in.data());
		transmit_type *recvptr = reinterpret_cast<transmit_type*>(out.data());
		const MPI_Datatype datatype = mpi::get_mpi_datatype<transmit_type>();
		int status = MPI_Allgatherv(sendptr, in_size, datatype, recvptr,
									sizes.data(), displacements.data(),
									datatype, comm);
		if (status != 0) 
		{
			std::cerr << "PE " << comm.rank() << ": MPI_Allgatherv returned "
			<< status << ", errno " << errno << std::endl;
			exit(-1);
		}
	}

	void GromacsHook::SyncToEngine()
	{
		gmxpush_();
	}

	void GromacsHook::SyncToSnapshot()
	{
		gmxpull_();
	}

	MPI_Comm GromacsHook::GetCommunicator()
	{
		return static_cast<MPI_Comm>(_snapshot->GetCommunicator());
	}

	template<typename T>
	void GromacsHook::PullToSSAGES(
		int iteration, 
		int natoms, 
		int* indices,
		int* types,
		T masses[], 
		T positions[][3], 
		T velocities[][3], 
		T forces[][3],
		T boxmat[3][3],
		double temperature,
		double potenergy,
		double kb)
	{
		// Resize vectors.
		auto& pos = _snapshot->GetPositions();
		pos.resize(natoms);
		auto& vel = _snapshot->GetVelocities();
		vel.resize(natoms);
		auto& frc = _snapshot->GetForces();
		frc.resize(natoms);
		auto& mass = _snapshot->GetMasses();
		mass.resize(natoms);
		auto& ids = _snapshot->GetAtomIDs();
		ids.resize(natoms);
		auto& typs = _snapshot->GetAtomTypes();
		typs.resize(natoms);

		// Reduce temperature/pressure/energy.
		auto& comm = _snapshot->GetCommunicator();

		// Atom weighted averages for temperature. Potential energy is extensive.
		// NOTE: pressure is just wrong.
		int ntot = 0;
		temperature *= natoms;
		MPI_Allreduce(&natoms, &ntot, 1, MPI_INT, MPI_SUM, comm);
		MPI_Allreduce(MPI_IN_PLACE, &temperature, 1, MPI_DOUBLE, MPI_SUM, comm);
		MPI_Allreduce(MPI_IN_PLACE, &potenergy, 1, MPI_DOUBLE, MPI_SUM, comm);
		temperature /= ntot;

		// Load em up.
		_snapshot->SetIteration(iteration);
		_snapshot->SetTemperature(temperature);
		_snapshot->SetEnergy(potenergy);
		_snapshot->SetKb(kb);

		for(int i = 0; i < natoms; ++i)
		{
			// Gromacs internally indexes atoms starting from 0.
			if(indices != nullptr)
				ids[i] = indices[i] + 1;
			else
				ids[i] = i + 1;
			typs[i] = types[i];

			mass[i] = masses[i];
			mass[i] = masses[i];
			mass[i] = masses[i];

			pos[i][0] = positions[i][0];
			pos[i][1] = positions[i][1];
			pos[i][2] = positions[i][2];

			vel[i][0] = velocities[i][0];
			vel[i][1] = velocities[i][1];
			vel[i][2] = velocities[i][2];

			frc[i][0] = forces[i][0];
			frc[i][1] = forces[i][1];
			frc[i][2] = forces[i][2];			
		}
		
		Matrix3 H;
		H << boxmat[0][0], boxmat[0][1], boxmat[0][2],
		     boxmat[1][0], boxmat[1][1], boxmat[1][2],
		     boxmat[2][0], boxmat[2][1], boxmat[2][2];
		_snapshot->SetHMatrix(H);

		std::vector<Vector3> gpos, gvel, gfrc;
		std::vector<double> gmass;
		std::vector<int> gids, gtyps; 
		
		allgatherv_serialize(comm, pos, gpos);
		allgatherv_serialize(comm, vel, gvel);
		allgatherv_serialize(comm, frc, gfrc);
		allgatherv_serialize(comm, mass, gmass);
		allgatherv_serialize<int,int>(comm, ids, gids);
		allgatherv_serialize<int,int>(comm, typs, gtyps);

		pos = gpos; 
		vel = gvel;
		frc = gfrc; 
		mass = gmass;
		ids = gids; 
		typs = gtyps;
	}

 	template<typename T>
	void GromacsHook::PushToGromacs(
		int natoms,
		int* indices, 
		int* types,
		T masses[],
		T positions[][3], 
		T velocities[][3],
		T forces[][3])
	{
		const auto& pos = _snapshot->GetPositions();
		const auto& vel = _snapshot->GetVelocities();
		const auto& frc = _snapshot->GetForces();
		const auto& mass = _snapshot->GetMasses();
		const auto& ids = _snapshot->GetAtomIDs();
		const auto& typs = _snapshot->GetAtomTypes();

		// Sync back only local atom data. all_gather guarantees 
		// sorting of data in vector based on rank. Compute offset
		// based on rank and sync back only the appropriate data.
		std::vector<int> allatoms;
		auto& comm = _snapshot->GetCommunicator();
		allatoms.reserve(comm.size());
		boost::mpi::all_gather(comm, natoms, allatoms);
		auto j = std::accumulate(allatoms.begin(), allatoms.begin() + comm.rank(), 0);
		
		for(int i = 0; i < natoms; ++i, ++j)
		{
			// Gromacs internally indexes atoms starting from 0.
			if(indices != nullptr)
				indices[i] = ids[j] - 1;
			
			types[i] = typs[j];

			masses[i] = mass[j];
			masses[i] = mass[j];
			masses[i] = mass[j];

			positions[i][0] = pos[j][0];
			positions[i][1] = pos[j][1];
			positions[i][2] = pos[j][2];

			velocities[i][0] = vel[j][0];
			velocities[i][1] = vel[j][1];
			velocities[i][2] = vel[j][2];

			forces[i][0] = frc[j][0];
			forces[i][1] = frc[j][1];
			forces[i][2] = frc[j][2];
		}
	}

// Explicit instantiation. Normally we'd let this happen in Gromacs, 
// but the boost mpi calls are giving lib dependency grief.
template void GromacsHook::PullToSSAGES<float>(
	int, 
	int, 
	int*,
	int*,
	float[],
	float[][3],
	float[][3],
	float[][3],
	float[3][3],
	double, 
	double,
	double);

template void GromacsHook::PushToGromacs<float>(
	int,
	int*, 
	int*, 
	float[],
	float[][3], 
	float[][3],
	float[][3]);

template void GromacsHook::PullToSSAGES<double>(
	int, 
	int, 
	int*,
	int*,
	double[],
	double[][3],
	double[][3],
	double[][3],
	double[3][3],
	double, 
	double,
	double);

template void GromacsHook::PushToGromacs<double>(
	int,
	int*, 
	int*, 
	double[],
	double[][3], 
	double[][3],
	double[][3]);
}
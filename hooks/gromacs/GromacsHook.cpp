#include "GromacsHook.h"
#include <boost/mpi.hpp>

namespace SSAGES
{
	void GromacsHook::SyncToEngine()
	{
		gmxpush_();
	}

	void GromacsHook::SyncToSnapshot()
	{
		gmxpull_();
	}

	template<typename T1, typename T2>
	void GromacsHook::PullToSSAGES(
		int iteration, 
		int natoms, 
		int* indices, 
		T1 masses[], 
		T2 positions[], 
		T2 velocities[], 
		T2 forces[],
		double temperature,
		double pressure,
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
		auto& types = _snapshot->GetAtomTypes();
		types.resize(natoms);

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
		_snapshot->SetPressure(pressure);
		_snapshot->SetEnergy(potenergy);
		_snapshot->SetKb(kb);
		for(int i = 0; i < natoms; ++i)
		{
			if(indices != nullptr)
				ids[i] = indices[i];
			else
				ids[i] = i;

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

		if(iteration % 1000 == 0 && comm.rank() == 0)
		{
			std::cout << "Iteration: " << iteration << " temperature " << temperature <<
			" pressure " << pressure << " potential " << potenergy << "\n"; 
		}
	}

 	template<typename T1, typename T2>
	void GromacsHook::PushToGromacs(
		int natoms,
		int* indices, 
		T1 masses[],
		T2 positions[], 
		T2 velocities[],
		T2 forces[])
	{
		const auto& pos = _snapshot->GetPositions();
		const auto& vel = _snapshot->GetVelocities();
		const auto& frc = _snapshot->GetForces();
		const auto& mass = _snapshot->GetMasses();
		const auto& ids = _snapshot->GetAtomIDs();

		for(int i = 0; i < natoms; ++i)
		{
			if(indices != nullptr)
				indices[i] = ids[i];

			masses[i] = mass[i];
			masses[i] = mass[i];
			masses[i] = mass[i];

			positions[i][0] = pos[i][0];
			positions[i][1] = pos[i][1];
			positions[i][2] = pos[i][2];

			velocities[i][0] = vel[i][0];
			velocities[i][1] = vel[i][1];
			velocities[i][2] = vel[i][2];

			forces[i][0] = frc[i][0];
			forces[i][1] = frc[i][1];
			forces[i][2] = frc[i][2];
		}
	}

// Explicit instantiation. Normally we'd let this happen in Gromacs, 
// but the boost mpi calls are giving lib dependency grief.
template void GromacsHook::PullToSSAGES<float, float[3]>(
	int, 
	int, 
	int*,
	float[],
	float[][3],
	float[][3],
	float[][3],
	double,
	double, 
	double,
	double);

template void GromacsHook::PushToGromacs<float, float[3]>(
	int,
	int*, 
	float[],
	float[][3], 
	float[][3],
	float[][3]);
}
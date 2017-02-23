#include "GromacsHook.h"
#include <boost/mpi.hpp>

using namespace boost;

namespace SSAGES
{
	void GromacsHook::SyncToEngine()
	{
		gmxpush_();
	}

	void GromacsHook::SyncToSnapshot()
	{
		gmxpull_();
		Hook::PostStepHook();
	}

	MPI_Comm GromacsHook::GetCommunicator()
	{
		return static_cast<MPI_Comm>(snapshot_->GetCommunicator());
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
		auto& pos = snapshot_->GetPositions();
		pos.resize(natoms);
		auto& vel = snapshot_->GetVelocities();
		vel.resize(natoms);
		auto& frc = snapshot_->GetForces();
		frc.resize(natoms);
		auto& mass = snapshot_->GetMasses();
		mass.resize(natoms);
		auto& ids = snapshot_->GetAtomIDs();
		ids.resize(natoms);
		auto& typs = snapshot_->GetAtomTypes();
		typs.resize(natoms);

		// Reduce temperature/pressure/energy.
		auto& comm = snapshot_->GetCommunicator();

		// Atom weighted averages for temperature. Potential energy is extensive.
		// NOTE: pressure is just wrong.
		int ntot = 0;
		temperature *= natoms;
		MPI_Allreduce(&natoms, &ntot, 1, MPI_INT, MPI_SUM, comm);
		MPI_Allreduce(MPI_IN_PLACE, &temperature, 1, MPI_DOUBLE, MPI_SUM, comm);
		MPI_Allreduce(MPI_IN_PLACE, &potenergy, 1, MPI_DOUBLE, MPI_SUM, comm);
		temperature /= ntot;

		// Load em up.
		snapshot_->SetIteration(iteration);
		snapshot_->SetTemperature(temperature);
		snapshot_->SetEnergy(potenergy);
		snapshot_->SetKb(kb);
		snapshot_->SetNumAtoms(natoms);

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
		snapshot_->SetHMatrix(H);
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
		const auto& pos = snapshot_->GetPositions();
		const auto& vel = snapshot_->GetVelocities();
		const auto& frc = snapshot_->GetForces();
		const auto& mass = snapshot_->GetMasses();
		const auto& ids = snapshot_->GetAtomIDs();
		const auto& typs = snapshot_->GetAtomTypes();

		for(int i = 0; i < natoms; ++i)
		{
			// Gromacs internally indexes atoms starting from 0.
			if(indices != nullptr)
				indices[i] = ids[i] - 1;
			
			types[i] = typs[i];

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
#pragma once 
#ifndef __GROMACS_HOOK_H
#define __GROMACS_HOOK_H

#include "Hook.h"
#include <functional>

namespace SSAGES
{
	class GromacsHook : public Hook
	{
	private:
		std::function<void()> gmxpush_, gmxpull_;

	protected: 
		// Implementation of the SyncToEngine interface.
		void SyncToEngine() override;

		// Implementation of the SyncToSnapshot interface. 
		void SyncToSnapshot() override;
	public: 
		GromacsHook() : gmxpush_(), gmxpull_()
		{} 

		// Get singleton instance of GromacsHook.
		static GromacsHook& Instance()
		{
			static GromacsHook instance;
			return instance;
		}

		void SetGMXPull(std::function<void()> func)
		{
			gmxpull_ = func;
		}

		void SetGMXPush(std::function<void()> func)
		{
			gmxpush_ = func;
		}

		void SyncToSSAGES()
		{
			SyncToSnapshot();
		}

		template<typename T1, typename T2>
		void PullToSSAGES(
			int iteration, 
			int natoms, 
			int* indices, 
			T1 masses[], 
			T2 positions[], 
			T2 velocities[], 
			T2 forces[])
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

			// Load em up.
			_snapshot->SetIteration(iteration);
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
 		}

 		template<typename T1, typename T2>
 		void PushToGromacs(
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

		~GromacsHook(){}
	};
}

#endif
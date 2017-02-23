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

		int niterations_;

		GromacsHook() : gmxpush_(), gmxpull_(), niterations_(0)
		{} 

	protected: 
		// Implementation of the SyncToEngine interface.
		void SyncToEngine() override;

		// Implementation of the SyncToSnapshot interface. 
		void SyncToSnapshot() override;
	public: 

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

		void SetPeriodicBoundaryConditions(int pbc)
		{
			// Definitions taken from Gromacs (enum.h).
			// enum {epbcXYZ, epbcNONE, epbcXY, epbcSCREW, epbcNR};
			switch(pbc)
			{
				case 0:
					snapshot_->SetPeriodicity({true, true, true});
					break;
				case 1:
					snapshot_->SetPeriodicity({false, false, false});
					break;
				case 2:
					snapshot_->SetPeriodicity({true, true, false});
					break;
				default:
					std::cerr << "Unsupported PBC specified in Gromacs." << std::endl;
					exit(-1);
			}
		}

		int GetIterationTarget() { return niterations_; }
		void SetIterationTarget(int niterations) { niterations_ = niterations; }

		MPI_Comm GetCommunicator();

		template<typename T>
		void PullToSSAGES(
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
			double kb);

 		template<typename T>
 		void PushToGromacs(
 			int natoms,
 			int* indices, 
 			int* types,
 			T masses[],
 			T positions[][3], 
 			T velocities[][3],
 			T forces[][3]);

		~GromacsHook(){}
	};
}

#endif
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
			T2 forces[],
			double temperature,
			double pressure,
			double potenergy,
			double kb);

 		template<typename T1, typename T2>
 		void PushToGromacs(
 			int natoms,
 			int* indices, 
 			T1 masses[],
 			T2 positions[], 
 			T2 velocities[],
 			T2 forces[]);

		~GromacsHook(){}
	};
}

#endif
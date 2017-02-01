#pragma once 

#ifndef OPENMD_VERSION_MAJOR
#include "config.h"
#include EXTRA_CONFIG
#endif

#include "Hook.h"
#include "brains/SimInfo.hpp"

namespace SSAGES
{
	class OpenMDHook : public Hook
	{
	private:
		OpenMD::SimInfo* siminfo_;

		// Private constructor to enforce singleton.
		OpenMDHook() {}
	
	protected: 
		// Implementation of the SyncToEngine interface.
		void SyncToEngine() override;

		// Implementation of the SyncToSnapshot interface. 
		void SyncToSnapshot() override;
	
	public: 

		// Get singleton instance of OpenMDHook.
		static OpenMDHook& Instance()
		{
			static OpenMDHook instance;
			return instance;
		}

		void SetSimInfo(OpenMD::SimInfo* siminfo)
		{
			siminfo_ = siminfo;
		}

		void SyncToSSAGES()
		{
			SyncToSnapshot();
		}
	};
}
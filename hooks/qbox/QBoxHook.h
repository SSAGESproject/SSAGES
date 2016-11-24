#pragma once 

#include "Hook.h"

namespace SSAGES
{
	class QBoxHook : public Hook
	{
	private: 
		int niterations_;
	protected:
		// Implementation of the SyncToEngine interface.
		void SyncToEngine() override
		{}

		// Implementation of the SyncToSnapshot interface. 
		void SyncToSnapshot() override
		{}
	public: 
		QBoxHook() 
		{}		
	};
}
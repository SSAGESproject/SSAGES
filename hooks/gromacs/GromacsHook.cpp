#include "GromacsHook.h"

// This cpp file is necessary to force gromacs into including necessary 
// SSAGES files.

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
}
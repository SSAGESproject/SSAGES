#ifdef FIX_CLASS

FixStyle(ssages,FixSSAGES)

#else
#ifndef LMP_FIX_SSAGES_H
#define LMP_FIX_SSAGES_H

#include "fix.h"
#include "Hook.h"

namespace LAMMPS_NS
{
	class FixSSAGES : public Fix, SSAGES::Hook
	{
	protected:
		void SyncToEngine() override;
		void SyncToSnapshot() override;

	public:
		FixSSAGES(class LAMMPS *, int, char**);
		void setup(int) override;
  		void pre_force(int) override;
  		void post_run() override;
  		int setmask() override;
	};
}

#endif
#endif
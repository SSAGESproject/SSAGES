#ifndef SSAGES_HOOMD_HOOK_H
#define SSAGES_HOOMD_HOOK_H

#include "Hook.h"
#include "hoomd/ExecutionConfiguration.h"
#include "hoomd/SystemDefinition.h"
#include "hoomd/HalfStepHook.h"
#include "hoomd/ComputeThermo.h"


// SSAGES Hook class for HOOMD
class SSAGESHOOMDHook : public SSAGES::Hook, public HalfStepHook
{
    private:
        // HOOMD Timestep
        unsigned int timestep_;

        // HOOMD SystemDefinition
        std::shared_ptr<SystemDefinition> sysdef_;

        // HOOMD ComputeThermo
        std::shared_ptr<ComputeThermo> thermo_;

    protected:
        // Implementation of the SyncToEngine interface.
        void SyncToEngine() override;

        // Implementation of the SyncToSnapshot interface.
        void SyncToSnapshot() override;

    public:
        SSAGESHOOMDHook();

        // Sets SystemDefinition to enable particle data access
        void setSystemDefinition(std::shared_ptr<SystemDefinition> sysdef)
        {
            sysdef_ = sysdef;
        }

        // Sets ComputeThermo to enable temperature and energy computation
        void setComputeThermo(std::shared_ptr<ComputeThermo> thermo)
        {
            thermo_ = thermo;
        }

        // Synchronize snapshot with SSAGES after computing forces
        void update(unsigned int timestep);

        ~SSAGESHOOMDHook();
};

#endif

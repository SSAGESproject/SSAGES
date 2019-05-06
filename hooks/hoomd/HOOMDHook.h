#ifndef SSAGES_HOOMD_HOOK_H
#define SSAGES_HOOMD_HOOK_H

#include "Hook.h"
#include "CVs/CVManager.h"

#include "hoomd/ExecutionConfiguration.h"
#include "hoomd/SystemDefinition.h"
#include "hoomd/HalfStepHook.h"
#include "hoomd/ComputeThermo.h"

#include "hoomd/extern/pybind/include/pybind11/pybind11.h"

// SSAGES Hook class for HOOMD
class HOOMDHook : public SSAGES::Hook, public HalfStepHook
{
    private:
        // HOOMD Timestep
        unsigned int timestep_;

        // HOOMD SystemDefinition
        std::shared_ptr<SystemDefinition> sysdef_;

        // HOOMD ComputeThermo
        std::shared_ptr<ComputeThermo> thermo_;

        // Python interpreter scope
        pybind11::object scope_;

    protected:
        // Implementation of the SyncToEngine interface.
        void SyncToEngine() override;

        // Implementation of the SyncToSnapshot interface.
        void SyncToSnapshot() override;

    public:
        HOOMDHook();

        //! Pre-simulation hook
        /*!
         * This should be called at the appropriate
         * time by the Hook implementation.
         */
        virtual void PreSimulationHook() override
        {
            // call base class method
            Hook::PreSimulationHook();

            // Set simulation engine hook on CVs
            for(auto& cv : cvmanager_->GetCVs())
            {
                cv->SetHook(*this);
            }
        }

        // Sets SystemDefinition to enable particle data access
        void setSystemDefinition(std::shared_ptr<SystemDefinition> sysdef) override
        {
            sysdef_ = sysdef;
        }

        // Get the SystemDefinition
        std::shared_ptr<SystemDefinition> getSystemDefinition()
            {
            return sysdef_;
            }

        // Sets ComputeThermo to enable temperature and energy computation
        void setComputeThermo(std::shared_ptr<ComputeThermo> thermo)
        {
            thermo_ = thermo;
        }

        // Synchronize snapshot with SSAGES after computing forces
        void update(unsigned int timestep) override;

        // Return the current time step
        unsigned int getTimeStep()
            {
            return timestep_;
            }

        // Set the python interpreter scope
        void setPyScope(pybind11::object scope)
            {
            scope_ = scope;
            }

        // Set the python interpreter scope
        pybind11::object getPyScope()
            {
            return scope_;
            }

        ~HOOMDHook();
};

#endif

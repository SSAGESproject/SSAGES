/**
 * This file is part of
 * SSAGES - Suite for Advanced Generalized Ensemble Simulations
 *
 * Copyright 2016 Yamil Colon <yamilcolon2015@u.northwestern.edu>
 *                Hythem Sidky <hsidky@nd.edu>
 *
 * SSAGES is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * SSAGES is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with SSAGES.  If not, see <http://www.gnu.org/licenses/>.
 */
#pragma once 

#include "CVs/CollectiveVariable.h"
#include "Validator/ObjectRequirement.h"
#include "Drivers/DriverException.h"
#include "Snapshot.h"
#include "schema.h"

#include "hoomd/ForceCompute.h"

namespace SSAGES
{
    //! Collective variable to calculate angle.
    class HOOMDWrapCV : public CollectiveVariable
    {
    private:
        //! The ForceCompute we wrap
        std::shared_ptr<ForceCompute> fc_;

        //! A reference to the HOOMD hook
        HOOMDHook& hook_;

    public:

        //! Constructor.
        /*!
         * \param fc The force compute to wrap
         *
         * Wrap a HOOMD ForceCompute into a CollectiveVariable
         *
         * \todo Bounds needs to be an input and periodic boundary conditions
         */
        HOOMDWrapCV(std::shared_ptr<ForceCompute> fc);
        fc_(fc)
        { }

        //! Initialize necessary variables.
        /*!
         * \param snapshot Current simulation snapshot.
         */
        void Initialize(const Snapshot& snapshot) override
        {
        #if 0
        using std::to_string;

        std::vector<int> found;
        snapshot.GetLocalIndices(atomids_, &found);
        int nfound = found.size();
        MPI_Allreduce(MPI_IN_PLACE, &nfound, 1, MPI_INT, MPI_SUM, snapshot.GetCommunicator());

        if(nfound != 3)
            throw BuildException({
                "TorsionalCV: Expected to find " + 
                to_string(3) + 
                " atoms, but only found " + 
                to_string(nfound) + "."
            });	
        #endif
        }

        //! This sets a reference to the simulation engine hook.
        /*!
         * This method is called during initialization and can be used by the CV
         * to access engine-specific data structures.
         */
        virtual void SetHook(const HOOMDHook& hook)
            {
            hook_ = hook;
            }

        //! Evaluate the CV.
        /*!
        * \param snapshot Current simulation snapshot.
        */
        void Evaluate(const Snapshot& snapshot) override
        {
            // evaluate the underlying ForceCompute
            fc_->compute(hook_.getTimeStep());

            // reduce CV
            val_ = fc->calcEnergySum();
        }

        //! Compute the particle forces and virial by scaling the gradient with a constant
        /*!
         * \param bias The constant scaling factor to apply to all particle forces
         * \param snapshot The snapshot that contains the forces
         *
         */
        virtual void ApplyBias(double bias, const Snapshot& snapshot) const
        {
            // evaluate the underlying ForceCompute (only if necessary)
            fc_->compute(hook_.getTimeStep());

            ArrayHandle<Scalar4> h_force(fc_->getForce(), access_location::host, access_mode::readwrite);
            ArrayHandle<Scalar> h_virial(fc_->getVirial(), access_location::host, access_mode::readwrite);
            ArrayHandle<Scalar4> h_torque(fc_->getTorque(), access_location::host, access_mode::readwrite);

            unsigned int virial_pitch = fc_->getVirial().getPitch();

            auto pdata = hook_->getSystemDefinition()->getParticleData();
            unsigned int n = pdata->getN();

            // multiply particle forces with bias
            for (unsigned int i = 0; i < n; ++i)
                h_force.data[i] *= bias;

            for (unsigned int j = 0; j < 6; ++j)
                for (unsigned int i = 0; i < n; ++i)
                    h_virial.data[j*virial_pitch + i] *= bias;

            #if 0
            // currently cannot set the external virial
            for (unsigned int j = 0; j < 6; ++j)
                fc_->setExternalVirial(j,fc_->getExternalVirial(j)*bias);
            #endif
        }

        //! \copydoc CollectiveVariable::BuildCV()
        static HOOMDWrapCV* Build(const Json::Value& json, const std::string& path)
        {
        Json::ObjectRequirement validator;
        Json::Value schema;
        Json::CharReaderBuilder rbuilder;
        Json::CharReader* reader = rbuilder.newCharReader();

        reader->parse(JsonSchema::HOOMDWrapCV.c_str(),
                      JsonSchema::HOOMDWrapCV.c_str() + JsonSchema::HOOMDWrapCV.size(),
                      &schema, NULL);
        validator.Parse(schema, path);

        // Validate inputs.
        validator.Validate(json, path);
        if(validator.HasErrors())
            throw BuildException(validator.GetErrors());

        std::vector<int> atomids;
        auto force = json["force"].asString();

        std::shared_ptr<ForceCompute> fc = context.attr(force).attr("cpp_force");
        return new HOOMDWrapCV(fc);
        }

        //! Returns true if this CV can modify the particle forces in the snapshot
        virtual bool ModifiesParticleForces() const
        {
        return false;
        }


    };
}

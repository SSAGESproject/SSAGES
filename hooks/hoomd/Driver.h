/**
 * This file is part of
 * SSAGES - Suite for Advanced Generalized Ensemble Simulations
 *
 * Copyright 2017 Hythem Sidky <hsidky@nd.edu>
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

#include <mpi.h>
#include "json/json.h"
#include "SSAGESHOOMDHook.h"
#include <hoomd/ExecutionConfiguration.h>


namespace py = pybind11;
using namespace py::literals;


namespace SSAGES
{

    //! Simulation driver.
    class Driver
    {
    private:
        //! Resource handler.
        class ResourceHandler* rh_;

		//! HOOMD hook.
		class std::shared_ptr<SSAGESHOOMDHook> hook_;

        //! HOOMD ExecutionConfiguration
        class std::shared_ptr<ExecutionConfiguration> exec_conf_;

        //! Number of steps to run in HOOMD
        unsigned int run_steps_;

        //! Command line arguments for hoomd.context.initialize()
        std::string cmd_args_;

    public:
        Driver(class ResourceHandler* rh, std::shared_ptr<SSAGESHOOMDHook> hook,
               std::shared_ptr<ExecutionConfiguration> exec_conf) :
        rh_(rh), hook_(hook), exec_conf_(exec_conf)
        {}

        //! Runs the Driver, which calls the python input file
        void Run();

        //! Build a new Driver from JSON.
		/*!
		 * \param json JSON root node containing driver (and children) specifications.
		 * \param world MPI communicator containing all processors.
		 * \return Pointer to newly created driver.
		 *
		 * \note Object lifetime is caller's responsibility!
		 */
		static Driver* Build(const Json::Value& json, const MPI_Comm& world);

        // Sets the number of steps to run
        void setRunSteps(unsigned int run_steps)
        {
            run_steps_ = run_steps;
        }

        // Sets the command line arguments for hoomd.context.initialize()
        void setCmdArgs(std::string cmd_args)
        {
            cmd_args_ = cmd_args;
        }

        ~Driver();
    };
}

/**
 * This file is part of
 * SSAGES - Software Suite for Advanced General Ensemble Simulations
 *
 * Copyright 2018 Bradley Dice <bdice@bradleydice.com>
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
#include <fstream>
#include <stdexcept>
#include "Driver.h"
#include "ResourceHandler.h"
#include <Python.h>
#include "hoomd/extern/pybind/include/pybind11/pybind11.h"
#include "hoomd/extern/pybind/include/pybind11/embed.h"
#include "hoomd/SystemDefinition.h"
#include "hoomd/Integrator.h"


namespace py = pybind11;
using namespace py::literals;
using namespace SSAGES;


namespace SSAGES
{

    void Driver::Run()
    {
        py::scoped_interpreter guard{}; // Start the interpreter and keep it alive
        py::module hoomd_module = py::module::import("hoomd");

        // Initialize HOOMD with the command line arguments
        hoomd_module.attr("context").attr("initialize")(cmd_args_);

        // Construct an execution context and set the HOOMD context exec_conf
        py::object exec_conf_obj = py::cast(exec_conf_);
        hoomd_module.attr("context").attr("exec_conf") = exec_conf_obj;

        // Evaluate the user script file
        py::eval_file(rh_->GetInput().c_str());
        py::object context = hoomd_module.attr("context").attr("current");

        // Get the current C++ SystemDefinition and update the HOOMDHook
        py::object sysdef_obj = context.attr("system_definition");
        std::shared_ptr<SystemDefinition> sysdef = sysdef_obj.cast<std::shared_ptr<SystemDefinition> >();
        hook_->setSystemDefinition(sysdef);

        // Get the current C++ Integrator and set the HOOMDHook
        py::object integrator_obj = context.attr("integrator").attr("cpp_integrator");
        std::shared_ptr<Integrator> integrator = integrator_obj.cast<std::shared_ptr<Integrator> >();
        integrator->setHalfStepHook(hook_);

        // Get the ComputeThermo for the entire system
        py::object group_all = hoomd_module.attr("group").attr("all")();
        py::object compute_thermo_obj = hoomd_module.attr("compute").attr("_get_unique_thermo")(group_all).attr("cpp_compute");
        hook_->setComputeThermo(compute_thermo_obj.cast<std::shared_ptr<ComputeThermo> >());

        // Add an Analyzer so that potential energy is always computed
        py::object log_func = hoomd_module.attr("analyze").attr("log");
        log_func(py::none(), py::make_tuple("potential_energy"), 1);

        // Run HOOMD
        hook_->update(0);  // this is a hack for initialization
        hook_->PreSimulationHook();
        hoomd_module.attr("run")(run_steps_);
        hook_->PostSimulationHook();
        return;
    }

    Driver* Driver::Build(const Json::Value& json, const MPI_Comm& world)
    {
        auto* rh = ResourceHandler::Build(json, world);


        std::string cmd_args;

        if (json.isMember("args"))
            {
            if (json["args"].isArray())
                {
                cmd_args = json["args"][0].asString();
                for (unsigned int i=1; i<json["args"].size(); i++)
                    {
                    cmd_args += " " + json["args"][i].asString();
                    }
                }
            else
                {
                cmd_args = json["args"].asString();
                }
            }
        else
            {
            cmd_args = "";
            }

        auto execution_mode = ExecutionConfiguration::executionMode::AUTO;

        if (cmd_args.find("--mode=cpu") != std::string::npos)
            {
            execution_mode = ExecutionConfiguration::executionMode::CPU;
            }
        else if (cmd_args.find("--mode=gpu") != std::string::npos)
            {
            execution_mode = ExecutionConfiguration::executionMode::GPU;
            }

        auto exec_conf = std::make_shared<ExecutionConfiguration>(
          execution_mode, std::vector<int>{-1}, false, false,
          std::shared_ptr<Messenger>(), 0, rh->GetLocalComm());

        auto hook = std::make_shared<HOOMDHook>();
        // This uses .get() on a shared_ptr because SSAGES doesn't accept
        // shared pointers. There may be a better way to do this
        rh->ConfigureHook(dynamic_cast<Hook*>(hook.get()));

        Driver* driver = new Driver(rh, hook, exec_conf);
        driver->setRunSteps(json["hoomd_steps"].asUInt());
        driver->setCmdArgs(cmd_args);
        return driver;
    }

    Driver::~Driver()
    {
        delete rh_;
    }
}

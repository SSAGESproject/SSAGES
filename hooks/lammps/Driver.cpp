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
#include "Driver.h"
#include "Hook.h"
#include "Drivers/DriverException.h"
#include "ResourceHandler.h"
#include "input.h"
#include "modify.h"
#include "fix.h"
#include <iostream>
#include <fstream>

using namespace LAMMPS_NS;

namespace SSAGES
{
    void Driver::Run()
    {
        lammps_ = new LAMMPS(0, NULL, rh_->GetLocalComm());

        // Execute file. 
        std::ifstream file(rh_->GetInput());
        std::string line;
        while(std::getline(file, line))
        {
            lammps_->input->one(line.c_str());
            if(line.find("ssages") != std::string::npos)
            {
                auto fid = lammps_->modify->find_fix("ssages");
                if(fid < 0)
                    throw BuildException({"Found SSAGES reference but could find fix!"});
                
                Hook* hook;
                if(!(hook = dynamic_cast<Hook*>(lammps_->modify->fix[fid])))
                    throw BuildException({"Unable to dynamic cast hook on node " + std::to_string(rh_->GetWalkerID())});
                rh_->ConfigureHook(hook);
            }
        }   
    }

    Driver* Driver::Build(const Json::Value& json, const MPI_Comm& world)
    {
        auto* rh = ResourceHandler::Build(json, world);
        return new Driver(rh);
    }
}
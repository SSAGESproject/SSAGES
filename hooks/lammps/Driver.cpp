/**
 * This file is part of
 * SSAGES - Software Suite for Advanced General Ensemble Simulations
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
#include "ResourceHandler.h"
#include "input.h"
#include "modify.h"
#include "fix.h"
#include <fstream>
#include <stdexcept>

using namespace LAMMPS_NS;

namespace SSAGES
{
    void Driver::Run()
    {
        lammps_ = new LAMMPS(0, NULL, rh_->GetLocalComm());

        Hook* hook;
        int fid = -1;
        std::string line;
        std::string linebuffer = "";
        size_t stridx;
        std::ifstream file(rh_->GetInput());
        
        // File doesn't exist or other error. 
        if(!file)
            throw std::runtime_error("Error opening file \"" + rh_->GetInput() + "\". Please check that the file exists.");
        
        // Execute file. 
        while(std::getline(file, line))
        {
            // If line(s) needs continuation, strip trailing
            // whitespace and "&" and stitch them together.
            stridx = line.find_last_not_of(" \t");
            if((stridx != std::string::npos) && (line.at(stridx) == '&'))
            {
                linebuffer.append(line,0,stridx);
            }
            else
            {
                line = linebuffer + " " + line;
                lammps_->input->one(line.c_str());
                
                // Only look for fix if it wasn't previously found.
                if(fid < 0)
                {
                    fid = lammps_->modify->find_fix("ssages");
                    if(fid >= 0)
                    {
                        if(!(hook = dynamic_cast<Hook*>(lammps_->modify->fix[fid])))
                            throw std::runtime_error("Unable to dynamic cast hook on node " + std::to_string(rh_->GetWalkerID()));
                        rh_->ConfigureHook(hook);
                    }
                }
                
                // Reset buffer
                linebuffer = "";
            }
        }
        
        // If we reach the end with no SSAGES fix then throw error.
        if(fid < 0) throw std::runtime_error("Could not find SSAGES fix.");
    }

    Driver* Driver::Build(const Json::Value& json, const MPI_Comm& world)
    {
        auto* rh = ResourceHandler::Build(json, world);
        return new Driver(rh);
    }

    Driver::~Driver()
    {
        delete rh_;
    }
}

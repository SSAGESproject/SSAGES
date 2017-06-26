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
#include "ResourceHandler.h"
#include "input.h"

using namespace LAMMPS_NS;

namespace SSAGES
{
    void Driver::Run()
    {
        lammps_ = new LAMMPS(0, NULL, rh_->GetLocalComm());
        lammps_->input->file(rh_->GetInput().c_str());
    }

    Driver* Driver::Build(const Json::Value& json, const MPI_Comm& world)
    {
        auto* rh = ResourceHandler::Build(json, world);
        return new Driver(rh);
    }
}
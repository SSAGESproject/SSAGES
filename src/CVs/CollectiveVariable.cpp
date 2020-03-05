/**
 * This file is part of
 * SSAGES - Software Suite for Advanced General Ensemble Simulations
 *
 * Copyright 2016 Ben Sikora <bsikora906@gmail.com>
 *                Emre Sevgen <sesevgen@uchicago.edu>
 *                Yamil Colon <yamilcolon2015@u.northwestern.edu>
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
#include "CollectiveVariable.h"
#include "AlphaRMSDCV.h"
#include "CVManager.h"
#include "AngleCV.h"
#include "AntiBetaRMSDCV.h"
#include "BoxVolumeCV.h"
#include "GyrationTensorCV.h"
#include "PairwiseCV.h"
#include "ParallelBetaRMSDCV.h"
#include "ParticleCoordinateCV.h"
#include "ParticlePositionCV.h"
#include "ParticleSeparationCV.h"
#include "RouseModeCV.h"
#include "TorsionalCV.h"
#include "RMSDCV.h"
#include "json/json.h"
#include <stdexcept>

namespace SSAGES
{
	CollectiveVariable* CollectiveVariable::BuildCV(const Json::Value &json, const std::string& path)
	{
		// Get CV type. 
		auto type = json.get("type", "none").asString();

		if(type == "Angle")
			return AngleCV::Build(json, path);
		else if(type == "BoxVolume")
			return BoxVolumeCV::Build(json, path);
		else if(type == "GyrationTensor")
			return GyrationTensorCV::Build(json, path);
		else if(type == "Pairwise")
			return PairwiseCV::Build(json, path);
		else if(type == "ParticleCoordinate")
			return ParticleCoordinateCV::Build(json, path);
		else if(type == "ParticlePosition")
			return ParticlePositionCV::Build(json, path);
		else if(type == "ParticleSeparation")
			return ParticleSeparationCV::Build(json, path);
		else if(type == "RouseMode")
			return RouseModeCV::Build(json, path);
		else if(type == "Torsional")
			return TorsionalCV::Build(json, path);
		else if (type == "AlphaRMSD")
			return AlphaRMSDCV::Build(json, path);
		else if (type == "ParallelBetaRMSD")
			return ParallelBetaRMSDCV::Build(json, path);
		else if (type == "AntiBetaRMSD")
			return AntiBetaRMSDCV::Build(json, path);
		else if (type == "RMSD")
			return RMSDCV::Build(json, path);
		else
			throw std::invalid_argument(path + ": Unknown CV type specified.");
	}

	std::map<std::string, unsigned int> CVManager::cvmap_ = {};
}

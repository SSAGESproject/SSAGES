#pragma once 

#include <iostream>

namespace SAPHRON
{
	class JsonSchema
	{
	public:
		//INSERT_DEC_HERE
		static std::string DirectorRestrictionC;
		static std::string Constraints;
		static std::string P2SAConnectivity;
		static std::string Connectivities;
		static std::string ParticleDistanceOP;
		static std::string WangLandauOP;
		static std::string ChargeFractionOP;
		static std::string RgOP;
		static std::string ElasticCoeffOP;
		static std::string Histogram;
		static std::string Simulation;
		static std::string DOSSimulation;
		static std::string LebwohlLasherFF;
		static std::string LennardJonesFF;
		static std::string DSFFF;
		static std::string LennardJonesTSFF;
		static std::string DebyeHuckelFF;
		static std::string ModLennardJonesTSFF;
		static std::string ForceFields;
		static std::string HardSphereFF;
		static std::string FENEFF;
		static std::string GayBerneFF;
		static std::string HarmonicFF;
		static std::string Worlds;
		static std::string SimpleWorld;
		static std::string Components;
		static std::string Particles;
		static std::string Site;
		static std::string Director;
		static std::string Selector;
		static std::string Blueprints;
		static std::string Position;
		static std::string Observers;
		static std::string DLMFileObserver;
		static std::string XYZObserver;
		static std::string JSONObserver;
		static std::string SpeciesSwapMove;
		static std::string AnnealChargeMove;
		static std::string VolumeSwapMove;
		static std::string Moves;
		static std::string AcidTitrationMove;
		static std::string DeleteParticleMove;
		static std::string InsertParticleMove;
		static std::string AcidReactionMove;
		static std::string RotateMove;
		static std::string FlipSpinMove;
		static std::string WidomInsertionMove;
		static std::string TranslateMove;
		static std::string TranslatePrimitiveMove;
		static std::string ParticleSwapMove;
		static std::string RandomIdentityMove;
		static std::string DirectorRotateMove;
		static std::string VolumeScaleMove;
		
	};
}
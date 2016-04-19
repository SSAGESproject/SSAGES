#include "Simulation.h"
#include "json/json.h"
#include "../Simulation/SimException.h"
#include "../Validator/ObjectRequirement.h"
#include "schema.h"
#include "SimInfo.h"
#include "StandardSimulation.h"
#include "DOSSimulation.h"

using namespace Json;

namespace SAPHRON
{
	Simulation* Simulation::BuildSimulation(const Value& json,
											WorldManager* wm, 
											ForceFieldManager* ffm,
											MoveManager* mm, 
											DOSOrderParameter* dop,
											Histogram* hist)
	{
		ObjectRequirement validator;
		Value schema; 
		Reader reader;

		Simulation* sim = nullptr;
	
		// Parse schema.
		reader.parse(JsonSchema::Simulation, schema);
		validator.Parse(schema, "#");

		// Validate input. 
		validator.Validate(json, "#");
		
		if(validator.HasErrors())
			throw BuildException(validator.GetErrors());

		// Set units. 
		auto units = json.get("units", "reduced").asString();
		if(units == "real")
		{
			auto& siminfo = SimInfo::Instance();
			siminfo.SetUnits(SimUnits::real);
		}

		// Sweeps.
		auto simtype = json["simtype"].asString();
		auto iterations = json["iterations"].asInt();
		auto mpi = json.get("mpi", 0).asInt();

		if(wm == nullptr || wm == NULL || 
		   ffm == nullptr || ffm == NULL ||
		   mm == nullptr || mm == NULL)
			throw BuildException({"#/simulation: All managers must not be null."});
		// Setup simulation. 
		if(simtype == "standard")
		{
			sim = new StandardSimulation(wm, ffm, mm);
		}
		else if(simtype == "DOS")
		{
			if(dop == nullptr || dop == NULL || 
			   hist == nullptr || hist == NULL)
				throw BuildException({"#/simulation: DOS simulation requires "
									  "initialized histogram and order parameter."});

			// Additional parameters for DOS.
			// Parse schema.
			reader.parse(JsonSchema::DOSSimulation, schema);
			validator.Parse(schema, "#");

			// Validate input. 
			validator.Validate(json, "#");

			auto resetfreq = json.get("reset_freq", 0).asInt();
			auto conv = json.get("convergence_factor", 1.0).asDouble();
			auto flatness = json.get("target_flatness", 0.80).asDouble();
			auto sync = json.get("sync_frequency", 100).asInt();
			auto equilib = json.get("equilibration", 0).asUInt();

			auto* dos = new DOSSimulation(wm, ffm, mm, dop, hist);
			dos->SetHistogramResetFrequency(resetfreq);
			dos->SetConvergenceFactor(conv);
			dos->SetTargetFlatness(flatness);
			dos->SetSyncFrequency(sync);
			dos->SetEquilibrationSweeps(equilib);

			sim = static_cast<Simulation*>(dos);
		}
		else
		{
			throw BuildException({"#/simulation/simtype: Unknown simtype \"" + simtype + "\"."});
		}

		sim->SetTargetIterations(iterations);
		if(mpi)
			sim->SetMovesPerIteration(mpi);

		return sim;
	}
}
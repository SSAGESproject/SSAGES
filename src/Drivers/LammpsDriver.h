#pragma once

#include "lammps.h"
#include "Drivers/Driver.h"
#include "../Validator/ObjectRequirement.h"
#include "../include/schema.h"

namespace mpi = boost::mpi;
using namespace LAMMPS_NS;
using namespace Json;

namespace SSAGES
{
	class LammpsDriver : public Driver 
	{
	private:

		//pointer to this local instance of lammps
		LAMMPS_NS::LAMMPS* _lammps;

		// The number of MD engine steps you would like to perform
		int _MDsteps;

		// The lammps logfile
		std::string _logfile;

	public:

		LammpsDriver(mpi::communicator& world_comm,
					 mpi::communicator& local_comm,
					 int walkerID,
					 Json::Value& jsonfile) : 
		Driver(world_comm, local_comm, walkerID, jsonfile), _lammps(), _MDsteps(), _logfile() 
		{
		};

		virtual void Run() override
		{
			std::string rline = "run " + std::to_string(_MDsteps);
			_lammps->input->one(rline.c_str());
		}

		// Run LAMMPS input file line by line and gather the fix/hook
		virtual void ExecuteInputFile(std::string contents) override
		{
			// Go through lammps.
			std::string token;
			std::istringstream ss(contents);
			while(std::getline(ss, token, '\n'))
				_lammps->input->one(token.c_str());

			auto fid = _lammps->modify->find_fix("ssages");
			if(!(auto* hook = dynamic_cast<Hook*>(_lammps->modify->fix[fid])))
			{
				if(_world.rank() == 0)
				{
					std::cerr << "Unable to dynamic cast hook. Error occurred" << std::endl;
					_world.abort(-1);			
				}
			}
		}

		virtual void BuildDriver() override
		{

			Value schema;
			ObjectRequirement validator;
			Reader reader;

			reader.parse(JsonSchema::LAMMPSDriver, schema);
			validator.Parse(schema, "#/Drivers");

			// Validate inputs.
			validator.Validate(_root, "#/Drivers");
			if(validator.HasErrors())
				throw BuildException(validator.GetErrors());

			_MDsteps = _root.get("MDSteps",1).asInt();

			// Silence of the lammps.
			char **largs = (char**) malloc(sizeof(char*) * 5);
			for(int i = 0; i < 5; ++i)
				largs[i] = (char*) malloc(sizeof(char) * 1024);
			sprintf(largs[0], " ");
			sprintf(largs[1], "-screen");
			sprintf(largs[2], "none");
			sprintf(largs[3], "-log");
			_logfile = _root.get("logfile", "none").asString();
			if(_logfile != "none")
				sprintf(largs[4], "%s-MPI_ID-%d",_logfile.c_str(), _wid);
			else
				sprintf(largs[4], "none");

			_lammps = std::make_shared<LAMMPS>(5, largs, MPI_Comm(_comm));

			// Free.
			for(int i = 0; i < 5; ++i)
				free(largs[i]);
			free(largs);

		}

		// Serialize
		virtual void Serialize(Json::Value& json) const override
		{
			json["MDSteps"] = _MDsteps;
			json["logfile"] = _logfile;
		}
	};
}
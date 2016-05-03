#pragma once

#include "lammps.h"
#include "Drivers/Driver.h"
#include "../Validator/ObjectRequirement.h"
#include "../include/schema.h"
#include "../Utility/BuildException.h"
#include "input.h"
#include "modify.h"
#include "fix.h"

namespace mpi = boost::mpi;
using namespace LAMMPS_NS;
using namespace Json;

namespace SSAGES
{
	class LammpsDriver : public Driver 
	{
	private:

		//pointer to this local instance of lammps
		std::shared_ptr<LAMMPS> _lammps;

		// The number of MD engine steps you would like to perform
		int _MDsteps;

		// The lammps logfile
		std::string _logfile;

	public:

		LammpsDriver(mpi::communicator& world_comm,
					 mpi::communicator& local_comm,
					 int walkerID) : 
		Driver(world_comm, local_comm, walkerID), _lammps(), _MDsteps(), _logfile() 
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

			std::cout<<token<<std::endl;

			auto fid = _lammps->modify->find_fix("ssages");
			if(fid < 0)
				throw BuildException({"Could not find ssages fix in given input file!"});

			if(!(_hook = dynamic_cast<Hook*>(_lammps->modify->fix[fid])))
			{
				throw BuildException({"Unable to dynamic cast hook on node " + std::to_string(_world.rank())});			
			}
		}

		virtual void BuildDriver(const Json::Value& json, const std::string& path) override
		{

			Value schema;
			ObjectRequirement validator;
			Reader reader;

			reader.parse(JsonSchema::LAMMPSDriver, schema);
			validator.Parse(schema, path);

			// Validate inputs.
			validator.Validate(json, path);
			if(validator.HasErrors())
				throw BuildException(validator.GetErrors());

			_MDsteps = json.get("MDSteps",1).asInt();
			_inputfile = json.get("inputfile","none").asString();
			
			// Silence of the lammps.
			char **largs = (char**) malloc(sizeof(char*) * 5);
			for(int i = 0; i < 5; ++i)
				largs[i] = (char*) malloc(sizeof(char) * 1024);
			sprintf(largs[0], " ");
			sprintf(largs[1], "-screen");
			sprintf(largs[2], "none");
			sprintf(largs[3], "-log");
			_logfile = json.get("logfile", "none").asString();
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
			json["type"] = "LAMMPS";
			json["number processors"] = _comm.size();
			if(_inputfile != "none")
				json["inputfile"] = _inputfile;

			// Need CVs and Methods still
			

		}
	};
}

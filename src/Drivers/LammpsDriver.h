#pragma once

#include "../JSON/Serializable.h"
#include "json/json.h"
#include <vector>
#include <boost/mpi.hpp>
#include "CVs/CollectiveVariable.h"
#include "Methods/Method.h"
#include "lammps.h"
#include <boost/mpi.hpp>

namespace mpi = boost::mpi;
using namespace LAMMPS_NS;

namespace SSAGES
{
	class LAMMPSDriver : public Driver 
	{
	private:

		auto _lammps = std::make_shared<LAMMPS>(5, largs, MPI_Comm(_comm));

		// The number of MD engine steps you would like to perform
		int _MDsteps;

		std::string _logfile;

	public:

		LammpsDriver(mpi::communicator& world_comm,
					 mpi::communicator& local_comm,
					 int walkerID,
					 Json::Value& jsonfile) : 
		Driver(world_comm, local_comm, walkerID, jsonfile) 
		{
		};

		virtual void Run() override
		{
			std::string rline = "run " + std::string(_MDsteps);
			lammps->input->one(rline.c_str());
		}

		// Run LAMMPS input file line by line and gather the fix/hook
		virtual void ExecuteInputFile(contents) override
		{
			// Go through lammps.
			std::string token;
			std::istringstream ss(contents);
			while(std::getline(ss, token, '\n'))
				lammps->input->one(token.c_str());

			// Get hook from lammps modify.
			// Horrid, I know.
			auto fid = lammps->modify->find_fix("ssages");
			if(!(auto* hook = dynamic_cast<Hook*>(lammps->modify->fix[fid])))
			{
				if(_comm.rank() == 0)
				{
					std::cerr << "Unable to dynamic cast hook on node "<<_wid<<". Error occurred" << std::endl;
					world.abort(-1);			
				}
			}
		}

		virtual void BuildDriver(const Json::Value& json, std::string path) override
		{

			Value schema;
			ObjectRequirement validator;
			Reader reader;

			reader.parse(JsonSchema::LAMMPSDriver, schema);
			validator.Parse(schema, path);

			// Validate inputs.
			validator.Validate(root, path);
			if(validator.HasErrors())
				throw BuildException(validator.GetErrors());

			_MDsteps = json.get("MDSteps",1).asInt();

			// Silence of the lammps.
			char **largs = (char**) malloc(sizeof(char*) * 5);
			for(int i = 0; i < 5; ++i)
				largs[i] = (char*) malloc(sizeof(char) * 1024);
			sprintf(largs[0], " ");
			sprintf(largs[1], "-screen");
			sprintf(largs[2], "none");
			sprintf(largs[3], "-log");
			_logfile = json.get("logfile", "none").asString();
			if(logfile != "none")
				sprintf(largs[4], "%s-MPI_ID-%d",logfile, _wid;
			else
				sprintf(largs[4], "none";

			_lammps = std::make_shared<LAMMPS>(5, largs, MPI_Comm(_comm));

		}

		// Serialize
		void Serialize(Json::Value& json) override
		{
			json["MDSteps"] = _MDsteps;
			json["logfile"] = _logfile;
		}
	};
}


			// Initialize snapshot. 
			Snapshot snapshot(walker, wid);

			hook->SetSnapshot(&snapshot);

			// Add methods and CV's here.
			///////Test Umbrella//////////////////////////////
			hook->AddListener(new Umbrella(world, walker, {std::stod(argv[3])}, {std::stod(argv[4])}, 1));
			hook->AddCV(new TorsionalCV(1, 5, 8, 11));
			// hook->AddCV(new ImproperCV(8, 5, 1, 11));

			//Mock method
			//hook->AddListener(new MockMethod(world, walker,1));

			///////Test MetaDynamics//////////////////////////
			//hook->AddListener(new Meta(0.5, {0.05, 0.05}, 500, 1));
			//hook->AddCV(new AtomCoordinateCV(1, 0));
			//hook->AddCV(new AtomCoordinateCV(1, 1));

			///////Test Elastic Band////////////////////////
			// Set up centers of each node based on toy system
			// if((int)world.size() < 3)
			//   {
			//     if(world.rank() == 0)
			//       std::cerr << "The elastic band method requires "
			// 		<< "at least 3 walkers." << std::endl;
			//     world.abort(-1);
			//   }

			// auto StartPointx = -1.1;
			// auto StartPointy = -1.05;
			// auto EndPointx = 1.1;
			// auto EndPointy = 1.15;

			// auto Nodediffxc = StartPointx + (int)world.rank()*(EndPointx - StartPointx)/(world.size()-1);
			// //		auto Nodediffyc = StartPointy + (EndPointy - StartPointy)*((double)world.rank()/(world.size()-1))*((double)world.rank()/(world.size()-1));
			// auto Nodediffyc = StartPointy + (EndPointy - StartPointy)*(int)world.rank()/(world.size()-1);

			// if((int)world.rank() + 1 == world.size())
			//   {
			//     Nodediffxc = EndPointx;
			//     Nodediffyc = EndPointy;
			//   }

			// hook->AddListener(new ElasticBand(world, walker,
			// 				  5000, 20000, 1000, 100, {Nodediffxc,Nodediffyc}, {100.0, 100.0}, 100.0, 0.001, 1));
			// hook->AddCV(new AtomCoordinateCV(1, 0));
			// hook->AddCV(new AtomCoordinateCV(1, 1));
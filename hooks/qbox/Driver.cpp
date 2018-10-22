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

#include <sys/stat.h>
#include <fstream>
#include <cstdio>
#include <unistd.h>
#include "json/json.h"
#include "Driver.h"
#include "QboxHook.h"
#include "ResourceHandler.h"

namespace SSAGES
{
	//! Waits for file to exist.
	//! \note Taken from twin.C in Qbox 1.63.5.
	void WaitForFile(std::string& lockfile)
	{
		struct stat statbuf; 
		int status; 
		do
		{
			// Status is 0 if file exists.
			status = stat(lockfile.c_str(), &statbuf);
			usleep(100000);
		}
		while( status != 0);
	}

	//! Waits for file to be deleted. 
	//! \note Taken from twin.C in Qbox 1.63.5.
	void WaitForNoFile(std::string& lockfile)
	{
		struct stat statbuf; 
		int status; 
		do
		{
			status = stat(lockfile.c_str(), &statbuf);
			usleep(100000);
		}
		while(status == 0);
	}

	void SendRunCommand(std::string& runfile, std::string& lockfile, int qmiter, int wfiter )
	{
		std::ofstream fin; 
		fin.open(runfile, std::ofstream::app);
		fin << "run 1 " << qmiter << " " << wfiter << std::endl;
		fin.close(); 
		remove(lockfile.c_str());
	}

	std::string Driver::CheckStorageFiles(std::string& storagefile,int run_no){ ////
		// This function check if there is already a backup file on disk
		// When it finds that the new backup file does not exist, it returns so it can be
		// open and user
		std::fstream file; ////
		file.open(storagefile, std::ios_base::out | std::ios_base::in); //// 
		std::string new_storage, return_storage; ////
		if(file.is_open()){////
			run_no++;////
			new_storage = "ssages_out_" + std::to_string(rh_->GetWalkerID()) + "_run_"+std::to_string(run_no)+".xml"; ////
			return_storage = CheckStorageFiles(new_storage,run_no);////
		}else{ ////
			return_storage = "ssages_out_" + std::to_string(rh_->GetWalkerID()) + "_run_"+std::to_string(run_no)+".xml"; ////
		}////
		file.close();////
		return return_storage;////
	}

	void Driver::Run()
	{
		using std::ofstream;
		using std::ifstream;

		// In and out files for Qbox communication. 
		auto infile = "ssages_in_" + std::to_string(rh_->GetWalkerID());
		auto outfile = "ssages_out_" + std::to_string(rh_->GetWalkerID());
		auto lockfile = infile + ".lock";

		std::cout << infile << std::endl ; 

		// Wait for initial lockfile to exist. This means Qbox is ready 
		WaitForFile(lockfile);
	
		// Write initial input file to qbox. 
		ofstream fin;
		fin.open(infile, std::ofstream::trunc);
		fin << rh_->GetInput() << std::endl;
		fin.close();

		// Remove lock file to get Qbox running.
		remove(lockfile.c_str());

		// Wait for Qbox to finish. 
		WaitForFile(lockfile);
		qbhook_->XMLToSSAGES(outfile);
		qbhook_->PreSimulationHook();

		// If first cycle save the first xml output ////
		int run_no = 0; ////
		std::string storage ; ////
		auto storage_zero = "ssages_out_" + std::to_string(rh_->GetWalkerID()) + "_run_"+std::to_string(run_no)+".xml"; ////
		storage = CheckStorageFiles(storage_zero,run_no); ////
		ofstream fstorage;////
		fstorage.open(storage);////
		fstorage.close();////
	
		//Initialize commands (defines "extforces" in qbox).
		qbhook_->InitializeCommands(infile);
		remove(lockfile.c_str());
		WaitForFile(lockfile);

		for(int i = 0; i < mditerations_; ++i)
		{
			qbhook_->SSAGESToCommands(infile);

			// Write run command to top it off and close file.
			SendRunCommand(infile, lockfile, qmiterations_, wfiterations_);

			// Wait for Qbox to finish. 
			WaitForFile(lockfile);
			qbhook_->XMLToSSAGES(outfile);
			qbhook_->PostIntegrationHook();

			// Append to storage file ////
			ofstream of_a(storage, std::ios_base::binary | std::ios_base::app); ////
			ifstream if_b(outfile, std::ios_base::binary); ////
			of_a.seekp(0, std::ios_base::end); ////
			of_a << if_b.rdbuf(); ////
		}
		
		qbhook_->SSAGESToCommands(infile);
		SendRunCommand(infile, lockfile, qmiterations_, wfiterations_);
		WaitForFile(lockfile);
		qbhook_->XMLToSSAGES(outfile);
		qbhook_->PostSimulationHook();	

		// Instructing Qbox that the run is finished ////
		// And make each walkers print a xml restart file ////
		fin.open(infile, std::ofstream::trunc); ////
		auto restart_xml = "restart_"+std::to_string(rh_->GetWalkerID())+".xml";
		fin << "save "<< restart_xml << std::endl ;
		fin << "quit" << std::endl;////
		fin.close();////
		remove(lockfile.c_str());////
	}

	Driver* Driver::Build(const Json::Value& json, const MPI_Comm& world)
	{
		auto* rh = ResourceHandler::Build(json, world);
		auto* hook = new QboxHook();
		rh->ConfigureHook(dynamic_cast<Hook*>(hook));

		auto iter = json["md_iterations"].asInt(); 
		auto qmiter = json.get("qm_iterations", 30).asInt();
		auto wfiter = json.get("wf_iterations", 0).asInt();

		return new Driver(rh, hook, iter, qmiter, wfiter);
	}   

	Driver::~Driver()
	{
		delete rh_;
		delete qbhook_;
	} 
}

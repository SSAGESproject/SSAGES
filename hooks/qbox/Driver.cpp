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

	void Driver::Run()
	{
		using std::ofstream;

		// In and out files for Qbox communication. 
		auto infile = "ssages_in" + std::to_string(0);
		auto outfile = "ssages_out" + std::to_string(0);
		auto lockfile = infile + ".lock";

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
		}
		
		qbhook_->SSAGESToCommands(infile);
		SendRunCommand(infile, lockfile, qmiterations_, wfiterations_);
		WaitForFile(lockfile);
		qbhook_->XMLToSSAGES(outfile);
		qbhook_->PostSimulationHook();	

		fin.open(infile, std::ofstream::trunc);
		fin << "quit" << std::endl;
		fin.close();
		remove(lockfile.c_str());
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

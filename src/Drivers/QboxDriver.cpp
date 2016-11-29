#include "QboxDriver.h"
#include <sys/stat.h>
#include <fstream>
#include <cstdio>

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

	void SendRunCommand(std::string& runfile, std::string& lockfile)
	{
		std::ofstream fin; 
		fin.open(runfile, std::ofstream::app);
		fin << "run 0 30" << std::endl;
		fin.close(); 
		remove(lockfile.c_str());
	}

	void QboxDriver::Run()
	{
		using std::ofstream;

		// In and out files for Qbox communication. 
		auto infile = "ssages_in" + std::to_string(_wid);
		auto outfile = "ssages_out" + std::to_string(_wid);
		auto lockfile = infile + ".lock";

		// Wait for initial lockfile to exist. This means Qbox is ready 
		WaitForFile(lockfile);
		
		// Write initial input file to qbox. 
		ofstream fin;
		fin.open(infile, std::ofstream::trunc);
		fin << _inputfile << std::endl;
		fin.close();

		// Remove lock file to get Qbox running.
		remove(lockfile.c_str());

		// Wait for Qbox to finish. 
		WaitForFile(lockfile);
		qbhook_->XMLToSSAGES(outfile);
		qbhook_->PreSimulationHook();
			
		for(int i = 0; i < mdsteps_; ++i)
		{
			qbhook_->SSAGESToCommands(infile);

			// Write run command to top it off and close file.
			SendRunCommand(infile, lockfile);

			// Wait for Qbox to finish. 
			WaitForFile(lockfile);
			qbhook_->XMLToSSAGES(outfile);
			qbhook_->PostIntegrationHook();
		}	

		qbhook_->SSAGESToCommands(infile);
		SendRunCommand(infile, lockfile);
		WaitForFile(lockfile);
		qbhook_->XMLToSSAGES(outfile);
		qbhook_->PostSimulationHook();
	}
}
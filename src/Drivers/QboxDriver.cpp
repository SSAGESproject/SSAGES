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

	void QboxDriver::Run()
	{
		using std::ofstream;

		// In and out files for Qbox communication. 
		auto infile = "ssages_in" + std::to_string(_wid);
		auto outfile = "ssages_out" + std::to_string(_wid);
		auto lockfile = infile + ".lock";

		// Wait for initial lockfile to exist. This means Qbox is ready 
		WaitForFile(lockfile);

		/*
		// Write to infile the input file name.
		ofstream fin;
		fin.open(infile.c_str(), std::ofstream::trunc);
		fin << _inputfile << std::endl;
		fin.close();

		// Remove lock file and wait for Qbox to finish.
		remove(lockfile.c_str());
		WaitForFile(lockfile);
		*/
		qbhook_->XMLToSSAGES(outfile);

		std::cout << "Qbox done!" << std::endl;


	}
}
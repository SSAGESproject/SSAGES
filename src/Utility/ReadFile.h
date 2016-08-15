#pragma once

#include <math.h>
#include <array>
#include <string>
#include <vector>
#include <iostream>
#include <fstream>

// Temp read file class. This class will be integrated/changed
// into a larger read file class that uses
// readxyz as a function. Can be expanded to include other
// file types.

namespace SSAGES
{
	class ReadFile
	{
	public:

		//! Constructor.
		ReadFile()
		{

		}

		//! Deconstructor
		~ReadFile()
		{

		}

		//! Read xyz file
		/*!
		 * \param FileName Name of xyz file to read in.
		 */
		static std::vector<std::array<double,4>> ReadXYZ(std::string FileName)
		{
			std::vector<std::array<double,4>> refcoord;
			int numatoms = 0;
			std::string comments = "";
			std::ifstream infile;
			infile.open(FileName, std::ios::in);
			if(!infile.is_open())
  				throw std::runtime_error("File " + FileName + " does not exist.");
			
			std::string ignore;
			
			std::getline(infile, ignore); // Get number of atoms

			numatoms = std::atoi(ignore.c_str());
			if(numatoms <= 0)
				throw std::runtime_error("Must be more than 0 atoms or invalid number atoms defined.");

			refcoord.resize(numatoms);

			std::getline(infile, comments); // Get comments

			int i = 0;
			std::string line;

			while (i < numatoms)
			{
			    std::getline(infile,line);
			    std::istringstream iss(line);
			    if (!(iss >> refcoord[i][0] >> refcoord[i][1] >> refcoord[i][2] >> refcoord[i][3])) 
			   		throw std::runtime_error("Bad file format for " + FileName + " for atom " + std::to_string(i+1)); 
				i++;
			}

			if(std::getline(infile, line) && !line.empty())
				throw std::runtime_error("Bad end line, possibly too many atoms defined.");

			if(i != numatoms)
				throw std::runtime_error("Number atoms specified does not match number of atoms defined");

			return refcoord;
		}
	};
}

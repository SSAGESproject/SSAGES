/**
 * This file is part of
 * SSAGES - Software Suite for Advanced General Ensemble Simulations
 *
 * Copyright 2016 Yamil Colon <yamilcolon2015@u.northwestern.edu>
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
#pragma once

#include <array>
#include <fstream>
#include <string>
#include <vector>

// Temp read file class. This class will be integrated/changed
// into a larger read file class that uses
// readxyz as a function. Can be expanded to include other
// file types.

namespace SSAGES
{
	//! Utility class to read file
	/*!
	 * Simple utility class to read in a trajectory file.
	 *
	 * Supported file types:
	 *
	 * * xyz
	 */
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
		 * \return Vector containing information stored in file.
		 *
		 * Read in a xyz file. The information will be returned as a vector of
		 * 4-element arrays. Each array corresponds to one atom and stores the
		 * following information: atom-ID, x-coordinate, y-coordinate,
		 * z-coordinate.
		 */
		static std::vector<std::array<double,4>> ReadXYZ(std::string FileName)
		{
			std::vector<std::array<double,4>> refcoord;
			int numatoms = 0;
			std::string comments = "";
			std::ifstream infile;
			infile.open(FileName, std::ios::in);
			if(!infile.is_open())
			{
				throw std::runtime_error("ReadFile.h: File " + FileName + " does not exist.");
			}

			std::string line;
			std::getline(infile, line); // Get number of atoms

			numatoms = std::stoi(line);
			if(numatoms <= 0)
			{
				throw std::runtime_error("Must be more than 0 atoms or invalid number atoms defined.");
			}

			refcoord.resize(numatoms);

			std::getline(infile, comments); // Get comments

			int i = 0;
			std::array<double,4> row;
			while(infile >> row[0] >> row[1] >> row[2] >> row[3])
			{
				if(i < numatoms) refcoord[i] = row;
				i++;
			}
			if(i != numatoms)
			{
				throw std::runtime_error(
					"ReadFile: Header specified " + std::to_string(numatoms) +
					" atoms, but read in " + std::to_string(i) + ". " +
					"Check the format of " + FileName + "."
				);
			}

			return refcoord;
		}
	};
}

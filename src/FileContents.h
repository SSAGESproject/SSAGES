/**
 * This file is part of
 * SSAGES - Suite for Advanced Generalized Ensemble Simulations
 *
 * Copyright 2016 Ben Sikora <bsikora906@gmail.com>
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

#include <string>
#include <cstdio>
#include <cerrno>

/*!
 * \file FileContents.h
 *
 * This file contains a collection of random helper functions.
 */
namespace SSAGES
{
	//! Read contents from a file
	/*!
	 * \param filename Name of the file to read from.
	 * \return String containing the contents of the file.
	 *
	 * Retrieves the contents of a file and returns them in a string. Throws
	 * exception on failure.
	 */
	inline std::string GetFileContents(const char *filename)
	{
		std::FILE *fp = std::fopen(filename, "rb");
		if (fp)
		{
			std::string contents;
			std::fseek(fp, 0, SEEK_END);
			contents.resize(std::ftell(fp));
			std::rewind(fp);

			// Stupid GCC bug. We do this to hide warnings.
			if(!std::fread(&contents[0], 1, contents.size(), fp))
				std::fclose(fp);
			else
				std::fclose(fp);

			return(contents);
		}
		throw(errno);
	}

	//! Gets file path from filename.
	/*!
	 * \param str String containing the absolute filename.
	 * \return String containing only the directory path.
	 */
	inline std::string GetFilePath(const std::string& str)
	{
		size_t found;
		found = str.find_last_of("/\\");
		if(found == str.npos)
			return "./";
		return str.substr(0, found);
	}
}
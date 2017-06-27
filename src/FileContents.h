/**
 * This file has been obtained from
 * SAPHRON - Statistical Applied PHysics through Random On-the-fly Numerics
 * https://github.com/hsidky/SAPHRON
 *
 * Copyright 2017 Hythem Sidky
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE. 
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
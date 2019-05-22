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

#include <fstream>
#include <string>

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
	inline std::string GetFileContents(const std::string filename)
	{
		std::ifstream ifs(filename);
		if(ifs) {
			std::string contents( (std::istreambuf_iterator<char>(ifs) ),
			                      (std::istreambuf_iterator<char>()    ) );
			return contents;
		}
		throw std::runtime_error("File " + filename + " does not exist.");
	}

	//! Gets file path from filename.
	/*!
	 * \param filename String containing the absolute filename.
	 * \return String containing only the directory path.
	 */
	inline std::string GetFilePath(const std::string filename)
	{
		size_t found;
		found = filename.find_last_of("/\\");
		if(found == filename.npos)
			return "./";
		return filename.substr(0, found);
	}
}
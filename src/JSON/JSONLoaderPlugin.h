/**
 * This file has been adapted from
 * SAPHRON - Statistical Applied PHysics through Random On-the-fly Numerics
 * https://github.com/hsidky/SAPHRON
 *
 * Copyright 2016 Hythem Sidky
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
#include <iostream>
#include <regex>
#include "FileContents.h"

namespace SSAGES
{
	//! Abstract class for JSON loader plugins.
	/*!
	 * \ingroup Core
	 */
	class JSONLoaderPlugin
	{
	public:
		//! Apply filter to string.
		/*!
		 * \param contents String whose contents will be modified.
		 * \param path Path for JSON path specification.
		 *
		 * This function will read in the contents string and apply custom
		 * filters to it.
		 */
		virtual void ApplyFilter(std::string& contents, const std::string& path) = 0;

		virtual ~JSONLoaderPlugin()
		{
		}
	};

	//! Class for JSON loader include plugin.
	/*!
	 * \ingroup Core
	 */
	class IncludePlugin : public JSONLoaderPlugin
	{
	public:
		//! Apply filter to string
		/*!
		 * \param contents String whose contents will be modified.
		 * \param path Path for JSON path specification.
		 *
		 * This filter replaces <tt>\@include(file.json)</tt> with
		 * contents.
		 */
		virtual void ApplyFilter(std::string& contents, const std::string& path) override
		{
			std::smatch matches;
			auto pattern = std::regex("\"@include\\((.*)\\)\"", std::regex::ECMAScript);
			while(regex_search(contents, matches, pattern))
			{
				for(size_t i = 1; i < matches.size(); ++i)
				{
					auto content = GetFileContents((path + "/" + matches[i].str()).c_str());
					auto rpattern = std::regex("\"@include\\(" +
											   matches[i].str() +
											   "\\)\"", std::regex::ECMAScript);
					contents = regex_replace(contents, rpattern, content);
				}
			}
		}

	};
}
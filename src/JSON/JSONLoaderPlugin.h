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
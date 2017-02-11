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
#include "json/json.h"
#include "JSONLoaderPlugin.h"
#include <fstream>
#include <memory>
#include <stdexcept>
#include <boost/mpi.hpp>

namespace mpi = boost::mpi;

namespace SSAGES
{
	//! Class for loading JSON content from files.
	/*!
	 * \ingroup Core
	 */
	class JSONLoader
	{
	private:
		//! List of plugins
		std::vector<std::unique_ptr<JSONLoaderPlugin>> plugins_;

	public:
		//! Constructor
		/*!
		 * The constructer automatically inserts the IncludePlugin into the
		 * list of plugins.
		 */
		JSONLoader() : plugins_(0)
		{
			plugins_.push_back(std::unique_ptr<JSONLoaderPlugin>(new IncludePlugin()));
		}

		//! Load file and return JSON tree.
		/*!
		 * \param filename Name of the JSON file to load.
		 * \param world MPI global communicator.
		 * \return JSON Value containing the contents of the file.
		 *
		 * This function loads JSON content from a given file. For each plugin
		 * in \c plugins_ the filter defined with this plugin will be applied
		 * to the contents loaded from the file.
		 */
		Json::Value LoadFile(const std::string& filename, boost::mpi::communicator& world)
		{
			std::string contents, path;
			if(world.rank() == 0)
			{
				contents = GetFileContents(filename.c_str());
				path = GetFilePath(filename);
			}

			mpi::broadcast(world, contents, 0);
			mpi::broadcast(world, path, 0);

			for(auto& plugin : plugins_)
				plugin->ApplyFilter(contents, path);

			// Read JSON.
			Json::Reader reader;
			Json::Value root;
			if(!reader.parse(contents, root))
				throw std::invalid_argument(
					filename + ": " + reader.getFormattedErrorMessages()
				);
			return root;
		}
	};
}
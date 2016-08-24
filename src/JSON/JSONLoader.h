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
		std::vector<std::unique_ptr<JSONLoaderPlugin>> _plugins;

	public:
		//! Constructor
		/*!
		 * The constructer automatically inserts the IncludePlugin into the
		 * list of plugins.
		 */
		JSONLoader() : _plugins(0)
		{
			_plugins.push_back(std::unique_ptr<JSONLoaderPlugin>(new IncludePlugin()));
		}

		//! Load file and return JSON tree.
		/*!
		 * \param filename Name of the JSON file to load.
		 * \param world MPI global communicator.
		 * \return JSON Value containing the contents of the file.
		 *
		 * This function loads JSON content from a given file. For each plugin
		 * in \c _plugins the filter defined with this plugin will be applied
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

			for(auto& plugin : _plugins)
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
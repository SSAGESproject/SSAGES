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
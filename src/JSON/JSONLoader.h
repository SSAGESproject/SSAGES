#pragma once 
#include "json/json.h"
#include "../Utils/Helpers.h"
#include "JSONLoaderPlugin.h"
#include <fstream>
#include <memory>
#include <stdexcept>

namespace SAPHRON
{

	// Class for loading JSON content from files.
	class JSONLoader
	{
	private:
		std::vector<std::unique_ptr<JSONLoaderPlugin>> _plugins;

	public:
		JSONLoader() : _plugins(0)
		{
			_plugins.push_back(std::unique_ptr<JSONLoaderPlugin>(new IncludePlugin()));
		}

		// Load file and return JSON tree.
		Json::Value LoadFile(const std::string& filename)
		{
			auto contents = GetFileContents(filename.c_str());
			auto path = GetFilePath(filename);

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
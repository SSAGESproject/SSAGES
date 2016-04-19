#pragma once 
#include "../Utils/Helpers.h"
#include <iostream>
#include <regex>

namespace SAPHRON
{
	// Abstract class for JSON loader plugins.
	class JSONLoaderPlugin
	{
	public:
		// Apply filter to string.
		virtual void ApplyFilter(std::string& contents, const std::string& path) = 0;
	};

	// Class for JSON loader include plugin. Replaces @include(file.json) with contents.
	class IncludePlugin : public JSONLoaderPlugin
	{
	public:
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
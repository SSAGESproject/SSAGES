#pragma once
#include <string>
#include <cstdio>
#include <cerrno>

// This a collection of random helper functions.
namespace SSAGES
{
	// Retrieves the contents of a file and returns them
	// in a string. Throws exception on failure.
	std::string GetFileContents(const char *filename)
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

	// Gets file path from filename.
	inline std::string GetFilePath(const std::string& str)
	{
		size_t found;
		found = str.find_last_of("/\\");
		if(found == str.npos)
			return "./";
		return str.substr(0, found);
	}
}
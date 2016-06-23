#pragma once

#include "../JSON/Serializable.h"
#include <vector>
#include <iostream>

// Forward declare.
namespace Json {
	class Value;
}

namespace SSAGES
{
	// Forward declare.
	class Grid;

	// Interface for a collective variable.
	class Grid
	{

	public:

		virtual ~Grid(){}

		// Return the nearest indices for given values.
		virtual std::vector<int> GetIndices(const std::vector<float> &values) const = 0;

		// Get the value at the current indices
		virtual float GetValue(const std::vector<int>& indices) const = 0;

		// Set the value at the current incices
		virtual void SetValue(const std::vector<int>& indices, float value) = 0;

		// Builds a Grid from a JSON node. Returns a pointer to the built grid.
		// If return value is nullptr, 
		// then an unknown error occurred. It will throw a BuildException on failure. 
		// Object lifetime is the caller's responsibility. 
		static Grid* BuildGrid(const Json::Value& json);

		// Overloaded function allowing JSON path specification.
		static Grid* BuildGrid(const Json::Value& json, 
							   const std::string& path);
	};
	
}
/**
 * This file is part of
 * SSAGES - Software Suite for Advanced General Ensemble Simulations
 *
 * Copyright 2017 Hythem Sidky <hsidky@nd.edu>
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

#include "CollectiveVariable.h"
#include "Drivers/DriverException.h"
#include "Validator/Requirement.h"
#include <vector>

namespace SSAGES
{
	//! Collective variable manager.
	/*!
	 * 
	 * CVManager is a class used to manage collective variables (CVs) and
	 * how they are exposed to methods. A metohd may wish to bias on a subset
	 * of available CVs. The CVManager provides a seamless interface to masking
	 * unwanted CVs and providing a suitable iterator which can be used to 
	 * iterate through the desired CVs. 
	 * 
	 * CVManager is also responsible for maintaining the lifetime of the CV 
	 * objects it contains. 
	 *
	 */
	class CVManager
	{
	private:
	//! List of collective variables.
	std::vector<CollectiveVariable*> cvs_;

	//! Map between CV names and ID's. 
	static std::map<std::string, unsigned int> cvmap_;

	public:
		CVManager() = default;

		//! Adds a CV to the CV manager.
		/*!
		 *
		 * \param cv Pointer to collective variable.
		 */
		void AddCV(CollectiveVariable* cv)
		{
			if(std::find(cvs_.begin(), cvs_.end(), cv) == cvs_.end())
				cvs_.push_back(cv);
		}

		//! Clears CVs from CV manager.
		/*
		 * \note This destroys all CVs stored in the CV manager!
		 */
		void ClearCVs()
		{
			for(auto& cv : cvs_)
				delete cv;
			cvs_.clear();
		}

		//! Get CV iterator.
		/*
		* \param mask Vector mask which contains the indices of
		*        which CV to include in the container.
		* \return Vector containing pointers to requested CVs.
		*/
		std::vector<CollectiveVariable*> GetCVs(const std::vector<unsigned int>& mask = std::vector<unsigned int>()) const
		{
			if(mask.empty())
				return cvs_;

			// Pack from mask.
			std::vector<CollectiveVariable*> cvs;
			for(auto& i : mask)
				cvs.push_back(cvs_[i]);

			return cvs;
		}

		//! Register CV name with map
		/*
		 * \param name Name of CV to register.
		 * \param id ID to associate with name.
		 *
		 * \note If a previous name is already used, it will override the old entry.
		 */
		static void AddCVtoMap(const std::string& name, unsigned int id)
		{
			cvmap_[name] = id;
		}

		//! Get CV id from map.
		/*
		 * \param cv JSON Value (CV) to look up.
		 * \param path Path for JSON path specification.
		 * \return ID of CV, -1 if nonexistent.
		 */
		static int LookupCV(const Json::Value& cv, const std::string& path)
		{
			int id = -1;
			if (cv.isString())
			{
				auto name = cv.asString();
				if(cvmap_.find(name) == cvmap_.end())
				{
					throw BuildException({path + ": CV mask name \"" + name + "\" does not exist."});
				}
				id = cvmap_.at(name);
			}
			else if(cv.isIntegral())
			{
				id = cv.asInt();
				if(id < 0 || id >= static_cast<int>(cvmap_.size()))
				{
					throw BuildException({path + ": CV mask index of " + std::to_string(id) + " does not exist. " +
					                      "Index must be nonnegative and less than " + std::to_string(cvmap_.size()) + "."});
				}
			}
			else
			{
				throw BuildException({path + ": CV mask must contain strings or nonnegative integers."});
			}

			return id;
		}

		~CVManager()
		{
			for(auto& cv : cvs_)
				delete cv;
		}
	};
}

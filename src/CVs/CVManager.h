/**
 * This file is part of
 * SSAGES - Suite for Advanced Generalized Ensemble Simulations
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

	public:
		CVManager() = default;

		//! Adds a CV to the CV manager. 
		/*!
		 * 
		 * \param cv Pointer to collective variable. 
		 */
		void AddCV(CollectiveVariable* cv)
		{
			cvs_.push_back(cv);
		}

		//! Get CV iterator. 
		 /* 
		  * \param mask Vector mask which contains the indices of
		  *        which CV to include in the container. 
		  * \return Vector containing pointers to requested CVs.
		  */
		std::vector<CollectiveVariable*> GetCVs(const std::vector<uint>& mask) const
		{
			if(mask.empty())
				return cvs_; 
			
			// Pack from mask. 
			std::vector<CollectiveVariable*> cvs;
			for(auto& i : mask)
				cvs.push_back(cvs_[i]);
			
			return cvs;
		}

		~CVManager()
		{
			for(auto& cv : cvs_)
            	delete cv;
		}
	};
}
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

#include "../Simulations/SimObserver.h"
#include "json/json.h"
#include <iostream>
#include <fstream>
#include <memory>

namespace SSAGES
{
	//! Observer with JSON output
	/*!
	 * JSON Observer class that creates a snapshot of the simulation and dumps
	 * it to a JSON file.
	 *
	 * \ingroup Core
	 */
	class JSONObserver : public SimObserver
	{
	private:
		std::unique_ptr<std::ofstream> _jsonfs; //!< Output file pointer.
		std::string _prefix; //!< Prefix for output file.
		Json::Value _root; //!< JSON root value.
		bool _writetoother; //!< Indicate which json to write to.

	public:
		//! Constructor
		/*!
		 * \param prefix Prefix for output file.
		 * \param world MPI global communicator.
		 * \param comm MPI local communicator.
		 * \param numwalkers Number of walkers.
		 * \param wid ID of walker to observe.
		 * \param frequency Frequency of observations.
		 */
		JSONObserver(const std::string& prefix, 
						boost::mpi::communicator world,
						boost::mpi::communicator comm,
						int numwalkers,
						int wid,
						unsigned int frequency = 1);

		//! Get Observer name.
		/*!
		 * \return String "JSON"
		 */
		virtual std::string GetName() const override{ return "JSON"; }

		//! Get filename prefix.
		/*!
		 * \return Filename prefix.
		 */
		std::string GetPrefix() const {return _prefix;}

		//! Visit driver
		/*!
		 * \param d Driver to visit.
		 */
		virtual void Visit(const Driver& d) override;

		//! \copydoc Serializable::Serialize()
		virtual void Serialize(Json::Value& json) const override;
		
		//! Operations before visiting.
		virtual void PreVisit() override;

		//! Operations after visition.
		virtual void PostVisit() override;

		//! Destructor.
		~JSONObserver()
		{
			if(_jsonfs)
				_jsonfs->close();
		}
	};
}
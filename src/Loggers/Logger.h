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

#include <mxx/comm.hpp>
#include <fstream>
#include "EventListener.h"

// Forward declare.
namespace Json {
	class Value;
}

namespace SSAGES
{
	//! Base class for logging SSAGES data. 
	/*!
	 * The base class for logging useful data in SSAGES 
	 * that is not necessarily written out by methods. Primarily 
	 * this includes the value(s) of the collective variable(s) 
	 * on the various walkers over time, the magnitude of the bias 
	 * if a method supports it, and so on. 
	 */
	class Logger : public EventListener
	{
	protected: 
		mxx::comm world_; //!< Global MPI communicator
		mxx::comm comm_; //!< Local MPI communicator

		//! Mask which identifies which CVs to log.
		std::vector<uint> cvmask_; 

		//! Name of logfile.
		std::string filename_; 

		//! Log file stream. 
		std::ofstream log_;

		//! Append mode? 
		bool append_;

	public: 
		//! Constructor
		/*!
		 * \param frequency Frequency of logging.
		 * \param filename File for logging.
		 * \param world Global MPI communicator.
		 * \param comm MPI communicator of walker.
		 * 
		 * \todo wid should be obtainable from somewhere. 
		 */
		Logger(uint frequency, const std::string& filename, const MPI_Comm& world, const MPI_Comm& comm) : 
		EventListener(frequency), world_(world), comm_(comm), cvmask_(), filename_(filename),
		append_(false)
		{}

		//! Logger call prior to simulation initiation.
		/*!
		 * \param snapshot Pointer to the simulation snapshot.
		 * \param cvmanager Collective variable manager.
		 *
		 * This function will be called before the simulation is started.
		 */
		virtual void PreSimulation(Snapshot* snapshot, const class CVManager& cvmanager) override;

		//! Logger call post integration.
		/*!
		 * \param snapshot Pointer to the simulation snapshot.
		 * \param cvmanager Collective variable manager.
		 *
		 * This function will be called after each integration step.
		 */
		virtual void PostIntegration(Snapshot* snapshot, const class CVManager& cvmanager) override;

		//! Logger call post simulation.
		/*!
		 * \param snapshot Pointer to the simulation snapshot.
		 * \param cvmanager Collective variable manager.
		 *
		 * This function will be called after the end of the simulation run.
		 */
		virtual void PostSimulation(Snapshot* snapshot, const class CVManager& cvmanager) override;
		
		//! Sets the collective variable mask.
		void SetCVMask(const std::vector<uint>& mask)
		{
			cvmask_ = mask;
		}

		//! Set append mode. 
		/*!
		 * \param append Whether to enable or disable append mode. 
		 */
		void SetAppend(bool append)
		{
			append_ = append;
		}

		//! Build a Logger from JSON node.
		/*!
		 * \param json JSON Value containing all input information.
		 * \param world MPI global communicator.
		 * \param comm MPI local communicator.
		 * \param path Path for JSON path specification.
		 * \return Pointer to the Logger built. nullptr if an unknown error occurred.
		 *
		 * \note Object lifetime is the caller's responsibility.
		 */
		static Logger* Build(const Json::Value& json, 
		                     const MPI_Comm& world, 
		                     const MPI_Comm& comm, 
		                     const std::string& path);

		//! Destructor
		virtual ~Logger() 
		{
		}	
	};
}
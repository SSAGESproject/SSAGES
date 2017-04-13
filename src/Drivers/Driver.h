/**
 * This file is part of
 * SSAGES - Suite for Advanced Generalized Ensemble Simulations
 *
 * Copyright 2016 Ben Sikora <bsikora906@gmail.com>
 *                Hythem Sidky <hsidky@nd.edu>
 *                Jiyuan Li <jyli@uchicago.edu>
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

#include <stdexcept>

#include "../JSON/Serializable.h"
#include "json/json.h"
#include <vector>
#include <boost/mpi.hpp>
#include <random>
#include "../Hook.h"
#include "../Methods/Method.h"
#include "../Constraints/Constraint.h"
#include "../Snapshot.h"
#include "../JSON/JSONLoader.h"
#include "../Simulations/SimObservable.h"
#include "../Simulations/SimObserver.h"
#include "../Observers/Visitable.h"
#include "../CVs/CollectiveVariable.h"

namespace SSAGES
{
	//! Abstract driver class for creating driver objects.
	/*!
	 * \ingroup Core
	 */
	class Driver : public SimObservable, public Serializable
	{

	protected:
		boost::mpi::communicator world_; //!< MPI global communicator
		boost::mpi::communicator comm_; //!< MPI local communicator

		//! The node id that this driver belongs to
		const int wid_;

		//! The hook that hooks into the MD simulation engine
		Hook* hook_;

		//! The snapshot of your system
		Snapshot* snapshot_;

		//! The Method that will be used
		Method* method_;

		//! The CVs that will be used
		CVList CVs_;

		//! Target number of iterations
		int iterations_;

		//! The observers that will be used for this driver
		ObserverList observers_;

		//! List of constraints to use
		ConstraintList constraints_;

		//! Local input file
		std::string inputfile_;

		//! MD engine Restart file name
		std::string restartname_;

		//! Read a restart file
		bool readrestart_;

		//! Random number generators for setting seeds
		std::random_device rd_;

		//! Alternative random number generator
		std::mt19937 gen_;

	public:
		//! Constructor
		/*!
		 * \param world MPI global communicator.
		 * \param comm MPI local communicator.
		 * \param walkerID The walker node for this driver.
		 */
		Driver(boost::mpi::communicator& world, 
			   boost::mpi::communicator& comm,
			   int walkerID) : 
		world_(world), comm_(comm), wid_(walkerID),
		hook_(nullptr), snapshot_(nullptr), method_(nullptr), CVs_(),
		iterations_(0), observers_(), inputfile_("none"), restartname_(),
		readrestart_(), rd_(), gen_(rd_())
		 {}

		//! Destructor
		virtual ~Driver()
		{
			for(auto& cv : CVs_)
				delete cv;

			for(auto& o : observers_)
				delete o;

			CVs_.clear();
			observers_.clear();

			delete snapshot_;

			delete method_;

			for(auto& constraint : constraints_)
				delete constraint;

			constraints_.clear();
		}

		//! Run simulation
		virtual void Run() = 0;

		//! Build the driver, which will create the hook and so forth
		/*!
		 * \param json JSON value containing input information.
		 * \param path Path for JSON path specification.
		 */
		virtual void BuildDriver(const Json::Value& json, const std::string& path) = 0;

		//! Read in driver specific input file (e.g. Lammps.in)
		/*!
		 * \param contents The contents of the input file.
		 *
		 * Execute the contents of an input file. The input file should be
		 * loaded with SSAGES::GetFileContents() first.
		 */
		virtual void ExecuteInputFile(std::string contents) = 0;

		//! Get the input file contents
		/*!
		 * \returns The contents of \c inputfile_.
		 *
		 * This function returns the contents of the input file.
		 *
		 * \note This function does not load the input file. If the input file
		 *       has not been loaded before, an empty string will be returned.
		 */
		std::string GetInputFile(){return inputfile_;}

		//! Set the input file
		/*!
		 * \param filename File name of the input file.
		 *
		 * This function sets the name of the input file for the driver.
		 */
		void SetInputFile(const std::string filename){ inputfile_ = filename;}

		//! Build CVs
		/*!
		 * \param json JSON Value containing input information.
		 * \param path Path for JSON path specification.
		 */
		void BuildCVs(const Json::Value& json, const std::string& path)
		{
			CollectiveVariable::BuildCV(json, CVs_, path);
		}

		// Build Constraints.
		/*!
		 * \param json JSON Value containing input information.
		 * \param path Path for JSON path specification.
		 */
		void BuildConstraints(const Json::Value& json, const std::string& path)
		{
			Constraint::BuildConstraint(json, constraints_, comm_, path);
		}

		//! Build method(s).
		/*!
		 * \param json JSON value containing input information.
		 * \param path Path for JSON path specification.
		 */
		void BuildMethod(const Json::Value& json, const std::string& path)
		{
			method_ = Method::BuildMethod(json, world_, comm_, path);
		}

		//! Build observer(s).
		/*!
		 * \param json JSON value containing input information.
		 * \param nwalks Number of walkers.
		 */
		void BuildObservers(const Json::Value& json,
			int nwalks)
		{
			SimObserver::BuildObservers(json, world_, comm_, nwalks, wid_, observers_);

			for(auto& o : observers_)
				this->AddObserver(o);
		}

		//! Finalize the setup
		/*!
		 * Create the snapshot and put all gathered values into the local hook.
		 * Set up listeners as well.
		 */
		void Finalize()
		{
			// Check that Hook is not a Null pointer
			if (!hook_) {
				throw std::runtime_error(
					"Trying to finalize simulation with invalid Hook."
				);
			}

			// Initialize snapshot. 
			snapshot_ = new Snapshot(comm_, wid_);
			snapshot_->SetTargetIterations(iterations_);

			/* Remove this later? */
			if(hook_ != nullptr)
			{
				// Set the hook to snapshot
				hook_->SetSnapshot(snapshot_);

				//Set the driver in the hook
				hook_->SetMDDriver(this);

			    if(method_)
				    hook_->AddListener(method_);

			    for(auto&c : constraints_)
				    hook_->AddListener(c);

				for(auto&cv : CVs_)
					hook_->AddCV(cv);
			}
		}

		//! \copydoc Serializable::Serialize()
		virtual void Serialize(Json::Value& json) const override
		{
			SerializeObservers(json["observers"]);

			method_->Serialize(json["method"]);

			auto& tmp = json["CVs"];
			for(unsigned int i = 0; i < CVs_.size();i++)
				CVs_[i]->Serialize(tmp[i]);

			json["inputfile"] = inputfile_;
			json["number processors"] = comm_.size();
			json["restart file"] = restartname_;
			json["read restart"] = readrestart_;
		}

     	// Accept a visitor.
		virtual void AcceptVisitor(Visitor& v) const override
		{
			v.Visit(*this);
		}
	};
}

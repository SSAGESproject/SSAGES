/**
 * This file is part of
 * SSAGES - Software Suite for Advanced General Ensemble Simulations
 *
 * Copyright 2016 Ashley Guo <azguo@uchicago.edu>
 *                Ben Sikora <bsikora906@gmail.com>
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

#include "StringMethod.h"
#include <fstream>

namespace SSAGES
{
	//! Finite Temperature Spring Method
	/*!
	 * \ingroup Methods
	 *
	 * Implementation of a multi-walker finite string
	 * method with hard wall voronoi cells and running block averages.
	 */
	class FiniteTempString : public StringMethod
	{
	private:	

		//! String modification parameter
		double kappa_;

		//! Number of steps to block average the CV's postions over
		unsigned int blockiterations_;

		//! Time step of string change
		double tau_;

		//! Minimum number of steps to apply umbrella sampling.
		unsigned int min_num_umbrella_steps_;

		//! Flag to run umbrella or not during post-integration
		int run_umbrella_;

		//! Iterator that keeps track of umbrella iterations
		unsigned int umbrella_iter_;

		//! Stores the last positions of the CVs
		std::vector<double> prev_CVs_;

		//! Checks if CV is in Voronoi cell
		/*!
		 * \param cvs List of CVs to check
		 *
		 * \return Boolean if all CVs are in the Voronoi cell
		 */
		bool InCell(const CVList& cvs) const;

		//! Updates the string according to the FTS method
		void StringUpdate() override;

		//! Flag for whether a system was to run umbrella sampling before checking against other systems
		bool reset_for_umbrella;

	public:
		//! Constructor
		/*!
		 * \param world MPI global communicator.
		 * \param comm MPI local communicator.
		 * \param centers List of centers.
		 * \param maxiterations Maximum number of iterations.
		 * \param blockiterations Number of iterations per block averaging.
		 * \param tau Value of tau (default: 0.1).
		 * \param cvspring Spring constants for cvs.
		 * \param kappa Value of kappa (default: 0.1).
		 * \param springiter Minimum number of umbrella steps.
		 * \param frequency Frequency with which this method is invoked.
		 *
		 * Constructs an instance of Finite String method.
		 */
		FiniteTempString(const MPI_Comm& world,
					     const MPI_Comm& comm,
					     const std::vector<double>& centers,
					     unsigned int maxiterations,
					     unsigned int blockiterations,
					     double tau,
					     const std::vector<double> cvspring,
					     double kappa,
					     unsigned int springiter,
			 		     unsigned int frequency) : 
		StringMethod(world, comm, centers, maxiterations, cvspring, frequency),
		kappa_(kappa), blockiterations_(blockiterations), tau_(tau), 
		min_num_umbrella_steps_(springiter), run_umbrella_(true),
		umbrella_iter_(1), reset_for_umbrella(false)
        {
			//! Store positions for starting trajectories
			prev_positions_.resize(1);

			//! Store velocities for starting trajectories
			prev_velocities_.resize(1);

			prev_IDs_.resize(1);

		}

		//! Post-integration hook.
		//! \copydoc Method::PostIntegration()
		void PostIntegration(Snapshot* snapshot, const class CVManager& cvmanager) override;

		//! Destructor
		~FiniteTempString() {}
	};
}
			

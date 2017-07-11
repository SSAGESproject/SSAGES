/**
 * This file is part of
 * SSAGES - Suite for Advanced Generalized Ensemble Simulations
 *
 * Copyright 2016 Hythem Sidky <hsidky@nd.edu>
 *                Jonathan K. Whitmer <jwhitme1@nd.edu>
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

#include "Method.h"
#include "Grids/Grid.h"
#include <fstream>
#include <vector>

namespace SSAGES
{
	//! Multidimensional hill
	/*!
	 * Structure representing a multidimensional hill (Gaussian) which is
	 * centered at "center" with widths "width" of height "height". A
	 * multidimensional Gaussian has one height but n centers and widths.
	 */
	struct Hill 
	{
		//! Hill center.
		std::vector<double> center;

		//! Hill width.
		std::vector<double> width;

		//! Hill height.
		double height;
		
		//! Constructs a multidimensional Hill (Gaussian)
		/*!
		 * \param center Hill center.
		 * \param sigma Hill width.
		 * \param height Hill height.
		 */
		Hill(const std::vector<double>& center, 
			 const std::vector<double>& sigma, 
			 double height) :
		 center(center), width(sigma), height(height)
		{}
	};

	//! "Vanilla" multi-dimensional Metadynamics.
	/*!
	 * Implementation of a "vanilla" multi-dimensional Metadynamics method with
	 * no bells and whistles.
	 *
	 * \ingroup Methods
	 */
	class Meta : public Method
	{
	private:	
		//! Hills.
		std::vector<Hill> hills_;

		//! Hill height.
		double height_;

		//! Hill widths.
		std::vector<double> widths_;

		//! Derivatives and temporary storage vectors.
		std::vector<double> derivatives_, tder_, dx_;

		//! Frequency of new hills
		unsigned int hillfreq_;

		//! CV Grid. 
		Grid<Vector>* grid_;

		//! Bounds 
		std::vector<double> upperb_, lowerb_; 

		//! Bound restraints. 
		std::vector<double> upperk_, lowerk_;

		//! Adds a new hill.
		/*!
		 * \param cvs List of CVs.
		 * \param iteration Current iteration.
		 */
		void AddHill(const CVList& cvs, int iteration);

		//! Computes the bias force.
		/*!
		 * \param cvs List of CVs.
		 */
		void CalcBiasForce(const CVList& cvs);

		//! Prints the new hill to file
		/*!
		 * \param hill Hill to be printed.
		 * \param iteration Current iteration.
		 */
		void PrintHill(const Hill& hill, int interation);

		//! Output stream for hill data.
		std::ofstream hillsout_;

	public:
		//! Constructor
		/*!
		 * \param world MPI global communicator.
		 * \param comm MPI local communicator.
		 * \param height Height of the hills to be deposited.
		 * \param widths Width of the hills to be deposited.
		 * \param hillfreq Frequency of depositing hills.
		 * \param frequency Frequency of invoking this method.
		 *
		 * Constructs an instance of Metadynamics method.
		 * \note The size of "widths" should be commensurate with the number of
		 *       CV's expected.
		 */
		Meta(const MPI_Comm& world,
			 const MPI_Comm& comm,
			 double height, 
			 const std::vector<double>& widths, 
			 const std::vector<double>& lowerb,
			 const std::vector<double>& upperb,
			 const std::vector<double>& lowerk,
			 const std::vector<double>& upperk,
			 Grid<Vector>* grid,
			 unsigned int hillfreq,
			 unsigned int frequency) : 
		Method(frequency, world, comm), hills_(), height_(height), widths_(widths), 
		derivatives_(0), tder_(0), dx_(0), hillfreq_(hillfreq), grid_(grid),
		upperb_(upperb), lowerb_(lowerb), upperk_(upperk), lowerk_(lowerk)
		{
		}

		//! Pre-simulation hook.
		/*!
		 * \param snapshot Current simulation snapshot.
		 * \param cvmanager Collective variable manager.
		 */
		void PreSimulation(Snapshot* snapshot, const class CVManager& cvmanager) override;

		//! Post-integration hook.
		/*!
		 * \param snapshot Current simulation snapshot.
		 * \param cvmanager Collective variable manager.
		 */
		void PostIntegration(Snapshot* snapshot, const class CVManager& cvmanager) override;

		//! Post-simulation hook.
		/*!
		 * \param snapshot Current simulation snapshot.
		 * \param cvmanager Collective variable manager.
		 */
		void PostSimulation(Snapshot* snapshot, const class CVManager& cvmanager) override;

		//! Load hills from file. 
		/*! 
		 * \param filename File name containing hills. 
		 *
		 * \note File format must match the output written by this method and 
		 *       the dimensionality of the hills must match the initialized 
		 *       Meta class. 
		 */
		void LoadHills(const std::string& filename);
		
		//! \copydoc Buildable::Build()
		static Meta* Build(const Json::Value& json, 
		                       const MPI_Comm& world,
		                       const MPI_Comm& comm,
					           const std::string& path);

		//! Destructor.
		~Meta() {}
	};
}
			

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
#include "../CVs/CollectiveVariable.h"
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
		std::vector<Hill> _hills;

		//! Hill height.
		double _height;

		//! Hill widths.
		std::vector<double> _widths;

		//! Gridding flag and grid
		bool _isgrid;

		//! Derivatives.
		std::vector<double> _derivatives;

		//! Bias magnitude.
		double _bias;

		//! Frequency of new hills
		unsigned int _hillfreq;

		//! Adds a new hill.
		/*!
		 * \param cvs List of CVs.
		 */
		void AddHill(const CVList& cvs);

		//! Computes the bias force.
		/*!
		 * \param cvs List of CVs.
		 */
		void CalcBiasForce(const CVList& cvs);

		//! Prints the new hill to file
		/*!
		 * \param hill Hill to be printed.
		 */
		void PrintHill(const Hill& hill);

		//! Output stream for hill data.
		std::ofstream _hillsout;

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
		Meta(boost::mpi::communicator& world,
			 boost::mpi::communicator& comm,
			 double height, 
			 const std::vector<double>& widths, 
			 unsigned int hillfreq,
			 bool isgrid,
			 unsigned int frequency) : 
		Method(frequency, world, comm), _hills(), _height(height), _widths(widths), 
		  _derivatives(0), _bias(0), _isgrid(isgrid), _hillfreq(hillfreq)
		{
		}

		//! Pre-simulation hook.
		/*!
		 * \param snapshot Current simulation snapshot.
		 * \param cvs List of CVs.
		 */
		void PreSimulation(Snapshot* snapshot, const CVList& cvs) override;

		//! Post-integration hook.
		/*!
		 * \param snapshot Current simulation snapshot.
		 * \param cvs List of CVs.
		 */
		void PostIntegration(Snapshot* snapshot, const CVList& cvs) override;

		//! Post-simulation hook.
		/*!
		 * \param snapshot Current simulation snapshot.
		 * \param cvs List of CVs.
		 */
		void PostSimulation(Snapshot* snapshot, const CVList& cvs) override;

		//! \copydoc Serializable::Serialize()
		/*!
		 * \warning Serialization not implemented yet!
		 */
		void Serialize(Json::Value& json) const override
		{

		}

		//! Destructor.
		~Meta() {}
	};
}
			

/**
 * This file is part of
 * SSAGES - Suite for Advanced Generalized Ensemble Simulations
 *
 * Copyright 2016 Ben Sikora <bsikora906@gmail.com>
 *                Joshua Lequieu <lequieu@uchicago.edu>
 *                Hadi Ramezani-Dakhel <ramezani@uchicago.edu>
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
#include "ForwardFlux.h"
#include "../CVs/CollectiveVariable.h"
#include <fstream>
#include <random>
#include <deque>
#include "../FileContents.h"
#include "../Drivers/DriverException.h"

namespace SSAGES
{
	//! ForwardFlux sampling method
	/*!
	 * \ingroup Methods
     * The notation used here is drawn largely from Allen, Valeriani and Rein ten Wolde. J. Phys.: Condens. Matter (2009) 21:463102. 
     * We recommend referring to this review if the reader is unfamiliar with the method, or our variable naming conventions.
	 */
	class DirectForwardFlux : public ForwardFlux
	{
	protected:

        //-----------------------------------------------------------------
        // Private Variables
        //-----------------------------------------------------------------

        //! Number of trials to attemt from each interface
        //! Note _M[0] sets the number of 'branches' for RBFFS and BGFFS?
        //std::vector<unsigned int> _M;

        //-----------------------------------------------------------------
        // Private Functions
        //-----------------------------------------------------------------
        

        //! Function that checks if interfaces have been crossed (different for each FFS flavor)
        void CheckForInterfaceCrossings(Snapshot*, const CVList&) override;

        //! Initialize the Queue
        void InitializeQueue(Snapshot*, const CVList&) override;

	public:
		//! Constructor
		/*!
		 * \param world MPI global communicator.
		 * \param comm MPI local communicator.
			 * \param frequency Frequency with which this method is invoked.
		 *
		 * Create instance of Forward Flux
		 */
		DirectForwardFlux(boost::mpi::communicator& world,
                    boost::mpi::communicator& comm, 
                    double ninterfaces, std::vector<double> interfaces,
                    unsigned int N0Target, std::vector<unsigned int> M,
                    bool initialFluxFlag, bool saveTrajectories, 
                    unsigned int currentInterface, unsigned int frequency)
		: ForwardFlux(world, comm, ninterfaces, interfaces, N0Target, M, 
					initialFluxFlag, saveTrajectories, currentInterface, frequency) {}

		//! Post-integration hook.
		/*!
		 * \param snapshot Current simulation snapshot.
		 * \param cvs List of CVs.
		 */
		void PostIntegration(Snapshot* snapshot, const CVList& cvs) override;

	};
}


/*
File Formats:
_indexfile
interface(some integer) dump_file_name(a string that contains interface and trial number)
example: 1 dump_1_10.xyz

dumpfile
atomid posx posy posz vx vy vz


*/

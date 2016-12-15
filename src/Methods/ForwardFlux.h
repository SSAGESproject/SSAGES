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
	class ForwardFlux : public Method
	{
	protected:

        //-----------------------------------------------------------------
        // Private Variables
        //-----------------------------------------------------------------

        //! Nested class to store different FFS Config IDs
        class FFSConfigID
        {
           public:
            unsigned int l; //!< Interface number
            unsigned int n;      //!< Configuration Number
            unsigned int a;      //!< Attempt number
            unsigned int lprev;      //!< Previous Interface number (i.e. traj I came from)
            unsigned int nprev;      //!< Previous Configuration Number
            unsigned int aprev;      //!< Previous Attempt number

            //! Constructor
            FFSConfigID(const unsigned int l, 
                        const unsigned int n, 
                        const unsigned int a, 
                        const unsigned int lprev, 
                        const unsigned int nprev, 
                        const unsigned int aprev): 
             l(l),n(n),a(a),lprev(lprev),nprev(nprev),aprev(aprev)
            {}

            ////! Yet another constructor
            //FFSConfigID(FFSConfigID in):
            //    l(in.l),n(in.n),a(in.a),lprev(in.lprev),nprev(in.nprev),aprev(in.aprev)
            //{}

            //! Empty Constructor
            FFSConfigID():
             l(0),n(0),a(0),lprev(0),nprev(0),aprev(0)
            {}

            ////! Overload = operator, not necessary, this is the default
            //FFSConfigID& operator=(const FFSConfigID& rhs){
            //    if (this == &rhs) return *this;
            //    else{
            //       l = rhs.l;
            //       n = rhs.n;
            //       a = rhs.a;
            //       lprev = rhs.lprev;
            //       nprev = rhs.nprev;
            //       aprev = rhs.aprev;
            //    }
            //}
        };

        //! random number generator
        std::default_random_engine _generator;

		//! Number of FFS interfaces
        //! note that _ninterfaces = n+1 where n is \lambda_n the interface defining B
		double _ninterfaces;

		//! FFS Interfaces
		std::vector<double> _interfaces;

        //! Current Interface
        unsigned int _current_interface;

		//! Previous cv position, used to determine if you've crossed an interface since last time
        double _cvvalue_previous;

		//!  current cv position
        double _cvvalue;

		//!  rate constant
        double _rate;

        //! Data structure that holds a Library N0 configurations at lambda0
        std::vector<FFSConfigID> Lambda0ConfigLibrary;

        //! Total Simulation Time spent in accumulating \ _N0
        double _N0TotalSimTime;
        
        //! Number of configurations to store at lambda0, target
        unsigned int _N0Target;

        //! Flux of trajectories out of state A. Denoted PhiA0 over h_A in Allen2009.
        double _fluxA0;

        //! Number of trials to attemt from each interface
        //! Note _M[0] sets the number of 'branches' for RBFFS and BGFFS?
        std::vector<unsigned int> _M;

        //! Number of attempts from interface i
        std::vector<unsigned int> _A;

        //! Flag to determine wheter fluxA0 should be calculated, seems not using this
        //bool _computefluxA0;

        //! Probability of going from lambda_{i} to lambda_{i+1}
        std::vector<double> _P;

        //! Number of successes from lambda_{i} to lambda_{i+1}
        //!  (might need to be 2d vector if multiple branches are used (with RBFFS)
        std::vector<unsigned int> _S;

		//! Current number of configurations currently stored at interface i
        //! This is somewhat redundant since _N[i] == _S[i-1], but for clarity 
		//! N[0] - current number of configurations collected at lambda0 (first interface) 
		std::vector<unsigned int> _N ;
        
        //! Keep track of jobs that have suceeded or failed but couldn't get reassigned a new task and must wait for the queue to get more jobs
        //! This could happen in DFFS once a job has finished but M[i] hasn't been reached (waiting on other jobs) 
        //! If this is the case I call it a 'zombie job', since the job is running, but isnt doing anything useful. Its just burning cpu cycles waiting for the queue to repopulate
        bool _pop_tried_but_empty_queue;

        //! if 1 compute initial flux
        bool _initialFluxFlag;
        bool initializeQueueFlag;

        //! The current FFSConfigID of this MPI process
        FFSConfigID myFFSConfigID;

        //! should the FFS trajectories be saved 
        bool _saveTrajectories;

        //! Counts the total number of failures
        //! eventually if I prune, will want this to be a vector, where it stores the number of failures at each interface
        //! however in the absence of pruning, a traj can only fail at lambda0, so this is just a scalar
        unsigned int _nfailure_total;

        //! commitor probability.
        //! The probability of a given configuration reaching B
        std::vector<std::vector<double>> _pB;


        /*! Queue
         *  When a given processor reaches an interface, it pulls a config from this Queue to figure out what it should do next
         *  This object should be syncronized between all FFS walkers (is walker the correct terminology here?)
         * technically this is a double-ended queue, this was mostly for debugging to allow element access of the queue (which std::queue doesn't allow). I use it like a queue though.
         */
        std::deque<FFSConfigID> FFSConfigIDQueue; 

        //! Directory of FFS output
        std::string _output_directory;

        //! file to which the current trajectory is written to
        std::ofstream _trajectory_file;

        //-----------------------------------------------------------------
        // Private Functions
        //-----------------------------------------------------------------
        
        //! Function that checks the initial structure that user provides.
        void CheckInitialStructure(const CVList&);

        //! Function to compute and write the initial flux
        void WriteInitialFlux();


        //! Function that adds new FFS configurations to the Queue
        //! Different FFS flavors can have differences in this method
        void AddNewIDsToQueue();

        //! Function checks if configuration has returned to A
        bool HasReturnedToA(double);

        //! Function checks if configuration has crossed interface specified since the last check
        /*! Simple function, given current and previous cv position, checks if interface i has been crossed. If crossed in positive direction, return +1, if crossed in negative direction return -1, if nothing crossed return 0
         */
        int HasCrossedInterface(double, double, unsigned int interface);

        //! Write a file corresponding to FFSConfigID from current snapshot
        void WriteFFSConfiguration(Snapshot *snapshot,FFSConfigID& ffsconfig, bool wassuccess);

        //! Read a file corresponding to a FFSConfigID into current snapshot
        void ReadFFSConfiguration(Snapshot *,FFSConfigID&,bool);
       
        //! Compute Initial Flux
        void ComputeInitialFlux(Snapshot*, const CVList&);

        //! Function that checks if interfaces have been crossed (different for each FFS flavor)
        virtual void CheckForInterfaceCrossings(Snapshot*, const CVList&) =0;

        //! Initialize the Queue
        virtual void InitializeQueue(Snapshot*, const CVList&) =0;

        //! Compute the probability of going from each lambda_i to lambda_{i+1} 
        /*!  
         *  Using number of successes and number of trials
         *  This will need to be different for each FFS flavor 
         */
        void ComputeTransitionProbabilities();

        //! Print the queue, useful for debugging
        void PrintQueue();

        //! Pop the queue, do MPI so that all procs maintain the same queue
        void PopQueueMPI(Snapshot*, const CVList&, bool);
        
        //! Compute the flux via brute force
        /*! Eventually this should be a new class that inherits from ForwardFlux, but for the time being I'll just hard code it
         *  This function takes the configurations at lambda0 and run them until they reach B (lambdaN) or return to A
         */
        void FluxBruteForce(Snapshot*, const CVList&);
        
        //! When simulation is finished, parse through the trajectories that reached B, and reconstruct the complete trajectory from where it started at A (lambda0)
        void ReconstructTrajectories(Snapshot *);

        //! When simulation is finished, recursively parse through the trajectories that reached B or failed back to A and calculate the Commitor Probability of that state going to B (_pB)
        void ComputeCommittorProbability(Snapshot *);
        
        //! Take the current config in snapshot and append it to the provided ofstream
        //! Current format is xyz style (including vx,vy,vz)
        void AppendTrajectoryFile(Snapshot*, std::ofstream&);

        //! Take the current config in snapshot and append it to the provided ofstream
        void OpenTrajectoryFile(std::ofstream&);

	public:
		//! Constructor
		/*!
		 * \param world MPI global communicator.
		 * \param comm MPI local communicator.
			 * \param frequency Frequency with which this method is invoked.
		 *
		 * Create instance of Forward Flux
		 */
		ForwardFlux(boost::mpi::communicator& world,
                    boost::mpi::communicator& comm, 
                    double ninterfaces, std::vector<double> interfaces,
                    unsigned int N0Target, std::vector<unsigned int> M,
                    bool initialFluxFlag, bool saveTrajectories,
                    unsigned int currentInterface, unsigned int frequency) : 
		 Method(frequency, world, comm), _ninterfaces(ninterfaces), _interfaces(interfaces), _N0Target(N0Target), 
         _M(M), _initialFluxFlag(initialFluxFlag), _saveTrajectories(saveTrajectories), _current_interface(currentInterface), _generator(1)
         {}

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
		virtual void PostIntegration(Snapshot* snapshot, const CVList& cvs) =0;

		//! Post-simulation hook.
		/*!
		 * \param snapshot Current simulation snapshot.
		 * \param cvs List of CVs.
		 */
		void PostSimulation(Snapshot* snapshot, const CVList& cvs) override;

		//! \copydoc Serializable::Serialize()
		void Serialize(Json::Value& json) const override
		{
			//Needed to run
			json["type"] = "ForwardFlux";

        }

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

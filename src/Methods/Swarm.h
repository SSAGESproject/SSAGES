/**        
- * This file is part of        
- * SSAGES - Suite for Advanced Generalized Ensemble Simulations        
- *     
- * Copyright 2016 Cody Bezik <bezik@uchicago.edu>      
- *     
- * SSAGES is free software: you can redistribute it and/or modify      
- * it under the terms of the GNU General Public License as published by        
- * the Free Software Foundation, either version 3 of the License, or       
- * (at your option) any later version.     
- *     
- * SSAGES is distributed in the hope that it will be useful,       
- * but WITHOUT ANY WARRANTY; without even the implied warranty of      
- * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the       
- * GNU General Public License for more details.        
- *     
- * You should have received a copy of the GNU General Public License       
- * along with SSAGES.  If not, see <http://www.gnu.org/licenses/>.     
- */
#include "StringMethod.h"
#include "../CVs/CollectiveVariable.h"
#include <fstream>
#include <iostream>

namespace SSAGES
{
    //! Swarm of Trajectories String Method
    /*!
     * \ingroup Methods
     *
     * Implementation of the swarm of trajectories string method.
     */
    class Swarm : public StringMethod //Calls Method constructor first, then Swarm constructor
    {
        private:

            //! Drift of CVs for one iteration
            std::vector<double> cv_drift_; 
                        
            //! Total number of MD steps for initialization for one iteration
            unsigned int initialize_steps_; 

            //! Length to run before harvesting a trajectory for unrestrained sampling
            unsigned int harvest_length_;

            //! Total number of restrained MD steps for one iteration
            unsigned int restrained_steps_; 

            //! Number of trajectories per swarm
            unsigned int number_trajectories_;

            //! Length of unrestrained trajectories
            unsigned int swarm_length_;

            //! Total number of unrestrained MD steps for one iteration
            unsigned int unrestrained_steps_;

            //! For indexing trajectory vectors
            unsigned int index_;

            //! Store atom IDs for starting trajecotires
            std::vector<Label> traj_atomids_; 

            //! Updates the positions of the string
            void StringUpdate();

            //! Helper function check if CVs are initialized correctly
            bool CVInitialized(const CVList& cvs);

            //! Flag for determing whether to perform initialization or not
            bool sampling_started;
            
            //! Flag for whether a snapshot was stored in the umbrella sampling
            bool snapshot_stored;        
            
            //! Flag for whether a system is initialized at a given iteration
            bool initialized;

            //! Flag for whether a system was initialized before it checked whether other systems were
            bool original_initialized; 

        public:

            //! Constructor.
            /*!
             * \param world MPI global communicator.
             * \param comm MPI local communicator.
             * \param centers List of centers.
             * \param maxiterations Maximum number of iterations.
             * \param cvspring Spring constants for CVs.
             * \param frequency Aplly method with this frequency.
             * \param InitialSteps Number of initial steps.
             * \param HarvestLength Length of trajectory before weighing.
             * \param NumberTrajectories Number of trajectories.
             * \param SwarmLength Lengt of the swarms.
             *
             * Constructs an instance of the swarm of trajectories method.
             */
            Swarm(boost::mpi::communicator& world, 
                    boost::mpi::communicator& comm, 
                    const std::vector<double>& centers,
                    unsigned int maxiterations,
                    const std::vector<double> cvspring,
                    unsigned int frequency,
                    unsigned int InitialSteps, 
                    unsigned int HarvestLength, 
                    unsigned int NumberTrajectories, 
                    unsigned int SwarmLength) : 
                StringMethod(world, comm, centers, maxiterations, cvspring, frequency),                 
                cv_drift_(), 
                initialize_steps_(InitialSteps), 
                harvest_length_(HarvestLength), 
                number_trajectories_(NumberTrajectories), 
                swarm_length_(SwarmLength)
        {
            cv_drift_.resize(centers_.size(), 0);
            prev_positions_.resize(number_trajectories_);
            prev_velocities_.resize(number_trajectories_);
            prev_IDs_.resize(number_trajectories_);
            //Additional initializing

            index_ = 0;  
            restrained_steps_ = harvest_length_*number_trajectories_; 
            unrestrained_steps_ = swarm_length_*number_trajectories_;
            sampling_started = false;
            snapshot_stored = false;

            iterator_ = 0; //Override default StringMethod.h initializing
        }

            //! Post-integration hook
            void PostIntegration(Snapshot* snapshot, const CVList& cvs) override; 

            void Serialize(Json::Value& json) const override
            {
                StringMethod::Serialize(json);

                json["flavor"] = "SWARM";

                json["initial_steps"] = initialize_steps_;
                json["harvest_length"] = harvest_length_;
                json["number_of_trajectories"] = number_trajectories_;
                json["swarm_length"] = swarm_length_; 
            }

            //! Destructor
            ~Swarm()
            {

            }
    };
}


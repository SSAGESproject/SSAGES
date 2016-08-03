#include "Method.h"
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
    class Swarm : public Method //Calls Method constructor first, then Swarm constructor
    {
        private:

            //! CV starting location values
            std::vector<double> _centers; 

            //! Drift of CVs for one iteration
            std::vector<double> _cv_drift; 

            //! Start of CV values for a trajectory
            std::vector<long double> _cv_start;
            
            //! Mesh index value for reparametrization
            double _alpha;  

            //! The computational node this image belongs to
            unsigned int _mpiid;

            //! Index of all CV values on the string
            std::vector<std::vector<double>> _worldstring;

            //! Total number of MD steps for initialization for one iteration
            unsigned int _initialize_steps; 

            //! Length to run before harvesting a trajectory for unrestrained sampling
            unsigned int _harvest_length;

            //! Total number of restrained MD steps for one iteration
            unsigned int _restrained_steps; 

            //! Number of trajectories per swarm
            unsigned int _number_trajectories;

            //! Length of unrestrained trajectories
            unsigned int _swarm_length;

            //! The local method iterator
            unsigned int _iterator;

            //! Total number of unrestrained MD steps for one iteration
            unsigned int _unrestrained_steps;

            //! Umbrella strength
            double _spring; 

            //! For indexing trajectory vectors
            int _index; 

            //! Output stream for string data
            std::ofstream _stringout;

            //! Store positions for starting trajectories
            std::vector<std::vector<Vector3>> _traj_positions;

            //! Store velocities for starting trajectories
            std::vector<std::vector<Vector3>> _traj_velocities;

            //! Number of nodes on a string
            unsigned int _numnodes;
            
            //! Updates the positions of the string
            void StringUpdate();

            //! Prints the string to a file
            void PrintString(const CVList& CV);

            //! Helper function to calculate distance
            double distance(std::vector<double>& x, std::vector<double>& y);

            //! Helper function check if CVs are initialized correctly
            bool CVInitialized(const CVList& cvs);

            //! Flag for determing whether to perform initialization or not
            bool sampling_started;

        public:

            //! Constructor
            /*!
             * Constructs an instance of the swarm of trajectories method.
             */
            Swarm(boost::mpi::communicator& world, 
                    boost::mpi::communicator& com, 
                    const std::vector<double>& centers, 
                    unsigned int NumNodes, 
                    double spring, 
                    unsigned int frequency, 
                    unsigned int InitialSteps, 
                    unsigned int HarvestLength, 
                    unsigned int NumberTrajectories, 
                    unsigned int SwarmLength) : 
                Method(frequency, world, com), 
                _centers(centers), 
                _cv_drift(), 
                _cv_start(), 
                _alpha(), 
                _mpiid(0), 
                _worldstring(), 
                _initialize_steps(InitialSteps), 
                _harvest_length(HarvestLength), 
                _number_trajectories(NumberTrajectories), 
                _swarm_length(SwarmLength), 
                _iterator(0), _spring(spring), 
                _numnodes(NumNodes)
        {
            
        }

            //! Pre-simulation hook
            void PreSimulation(Snapshot* snapshot, const CVList& cvs) override;

            //! Post-integration hook
            void PostIntegration(Snapshot* snapshot, const CVList& cvs) override; 

            //! Post-simulation hook
            void PostSimulation(Snapshot* snapshot, const CVList& cvs) override; 

            //! Set the method iteration.
            /*!
             * \param iter New value of the iteration.
             */
            void SetIteration(const int iter)
            {
                _iteration = iter;
            }		

            void Serialize(Json::Value& json) const override
            {
                json["type"] = "Swarm";
                for(size_t i = 0; i < _centers.size(); i++)
                {
                    json["centers"].append(_centers[i]);
                }

                json["number of nodes"] = _numnodes;
                json["spring"] = _spring;
                json["initial steps"] = _initialize_steps;
                json["harvest length"] = _harvest_length;
                json["number of trajectories"] = _number_trajectories;
                json["swarm length"] = _swarm_length;

                json["iteration"] = _iteration; 
            }

            //! Destructor
            ~Swarm()
            {

            }
    };
}


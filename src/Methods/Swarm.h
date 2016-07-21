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

            //! Number of iterations run so far
            unsigned int _currentiter;

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

            //! Store forces for starting trajectories
            std::vector<std::vector<Vector3>> _traj_forces;

            //! Number of nodes on a string
            unsigned int _numnodes;

            //! Debugging - drift scale parameter
            double _drift_scale; 
            
            //! Updates the positions of the string
            void StringUpdate();

            //! Prints the string to a file
            void PrintString(const CVList& CV);

        public:

            //! Constructor
            /*!
             * Constructs an instance of the swarm of trajectories method.
             */
            Swarm(boost::mpi::communicator& world, boost::mpi::communicator& com, const std::vector<double>& centers, unsigned int NumNodes, double spring, unsigned int frequency, unsigned int InitialSteps, unsigned int HarvestLength, unsigned int NumberTrajectories, unsigned int SwarmLength) : Method(frequency, world, com), _centers(centers), _cv_drift(), _cv_start(), _alpha(), _mpiid(0), _worldstring(), _currentiter(0), _initialize_steps(InitialSteps), _harvest_length(HarvestLength), _number_trajectories(NumberTrajectories), _swarm_length(SwarmLength), _iterator(0), _spring(spring), _numnodes(NumNodes), _drift_scale(1)
        {
            /*for(size_t i = 0; i < centers.size(); i++)
            {
                std::cout << _centers[i] << " " << centers[i] << std::endl;
            }*/
        }

            //! Pre-simulation hook
            void PreSimulation(Snapshot* snapshot, const CVList& cvs) override;

            //! Post-integration hook
            void PostIntegration(Snapshot* snapshot, const CVList& cvs) override; 

            //! Post-simulation hook
            void PostSimulation(Snapshot* snapshot, const CVList& cvs) override; 

            void Serialize(Json::Value& json) const override
            {

            }

            //! Destructor
            ~Swarm()
            {

            }
    };
}


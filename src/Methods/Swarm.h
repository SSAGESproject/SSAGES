#include "Method.h"
#include "../CVs/CollectiveVariable.h"

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

            //! Number of iterations run so far
            unsigned int _currentiter;

            //! The computational node this image belongs to
            unsigned int _mpiid;

            //! Number of nodes on a string
            unsigned int _numnodes;

            //! Length of initialization
            unsigned int _initialize_steps; 

            //! Length of restrained sampling 
            unsigned int _restrained_steps; 

            //! Number of trajectories per swarm
            unsigned int _number_trajectories;

            //! Length of unrestrained trajectories
            unsigned int _unrestrained_steps; 

            //! The local method iterator
            unsigned int _iterator;

            //! CV starting location values
            std::vector<double> _centers; 

            //! Drift of CVs for one iteration
            std::vector<double> _cv_drift; 

            //! Start of CV values for a trajectory
            std::vector<double> _cv_start;

            //! Mesh index value for reparametrization
            double _alpha; 

            //! Index of all CV values on the string
            std::vector<std::vector<double>> _worldstring;

            //! Umbrella strength
            double _spring; 

            //! For indexing trajectory vectors
            int _index; 

            //! Output stream for string data
            std::ofstream _stringout; 

            //! Updates the positions of the string
            void StringUpdate();

            //! Prints the string to a file
            void PrintString(const CVList& CV);

        public:

            //! Constructor
            /*!
             * Constructs an instance of the swarm of trajectories method.
             */
            Swarm(boost::mpi::communicator& world, boost::mpi::communicator& com, const std::vector<double>& centers, unsigned int NumNodes, double spring, unsigned int frequency, unsigned int InitialSteps, unsigned int RestrainedSteps, unsigned int NumberTrajectories, unsigned int UnrestrainedSteps) : Method(frequency, world, com) _centers(centers), _cv_prev(), alpha(), _mpiid(0), _worldstring(), _currentiter(0), _initiaize_steps(InitialSteps), _restrained_steps(RestrainedSteps), _number_trajectories(NumberTrajectories), _unrestrained_steps(UnrestrainedSteps), _iterator(0); 

            //! Pre-simulation hook
            void PreSimulation(Snapshot* snapshot, const CVList& cvs) override;

            //! Post-integration hook
            void PostIntegration(Snapshot* snapshot, const CVList& cvs) override; 

            //! Post-simulation hook
            void PostSimulation(Snapshot* snapshot, const CVList& cvs) override; 

            void Serialize(Json::Value& json) const override; 
            {

            }

            //! Destructor
            ~Swarm()
            {

            }
    };
}


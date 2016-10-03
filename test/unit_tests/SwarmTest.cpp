#include "gtest/gtest.h"
#include <boost/mpi.hpp>

#include "../src/Snapshot.h"


#define private public
#define protected public 

#include "../src/Methods/Swarm.h"
#include "../src/CVs/MockCV.h"

using namespace SSAGES;

// Test calculation up to accuracy of 10^-5
const double eps = 0.00001;

class SwarmTest : public ::testing::Test 
{
protected:

    virtual void SetUp() 
	{

		// Set up dummy CVs.
		CV1 = new MockCV(-0.5,{1,0,0},-2,2); // Initialized CV
		CV2 = new MockCV(-0.5,{0,1,0},-2,2); // Initialized CV
        CV3 = new MockCV(1.5,{1,0,0},-2,2); // Uninitialized CV

		cvlist_initialized.push_back(CV1);
		cvlist_initialized.push_back(CV2);
        cvlist_not_initialized.push_back(CV1);
        cvlist_not_initialized.push_back(CV3);
        
        mpiid = world.rank();

        // Set up Swarm parameters.		
        std::vector<double> _centers;
        _centers.resize(2); //Centers will contain two elements to check against
		if(mpiid == 0)
			_centers = {-0.49, -0.51};
		else if(mpiid == 1)
			_centers = {2.2, -1.8};
		else if(mpiid == 2)
			_centers = {2.8, 2.2};
        std::vector<double> _cvspring;
        
        unsigned int iteration = 5;
		unsigned int numnodes = world.size();
		worldstring = {{-0.49, -0.51}, {2.2, -1.8}, {2.8, 2.2}};

        //Dummy Swarm constructor - world, comm, centers, maxiterations, cvspring, frequency, initial steps, harvest length, number trajectories, swarm length
        Swarm_Method = new Swarm(world, comm, _centers, 1000, _cvspring, 1, 10000, 10, 100, 20); 
        Swarm_Method->_worldstring = worldstring;
        Swarm_Method->_mpiid = mpiid;
        Swarm_Method->_numnodes = numnodes;
        Swarm_Method->_iteration = iteration;
    }

	virtual void TearDown() 
	{
		delete CV1;
		delete CV2;
		delete CV3;

        delete Swarm_Method;
	}

    boost::mpi::communicator world;
    boost::mpi::communicator comm = world.split(world.rank() < 3 ? world.rank() : 4);
    
    unsigned int mpiid;
    std::vector<std::vector<double>> worldstring;

	MockCV* CV1;
	MockCV* CV2;
	MockCV* CV3;

	CVList cvlist_initialized;
    CVList cvlist_not_initialized;
    
    Swarm* Swarm_Method;
};

TEST_F(SwarmTest,dummytest)
{
    //The intention here is to suppress output from functions
    std::cout.setstate(std::ios_base::failbit); 

	// only run this test with 3 ~faux walkers~
	ASSERT_TRUE(3 == world.size()) << "ERROR: Swarm unit test must be run with 3 processes";

    //std::cout << " Hello 1 " << std::endl; //Debugging
    //Test for CV Initialized
    if(mpiid == 0)
    {
        bool _should_be_initialized = Swarm_Method->CVInitialized(cvlist_initialized);
        bool _should_not_be_initialized = Swarm_Method->CVInitialized(cvlist_not_initialized);

        //Backwards logic could be cleaned up in swarm itself
        //std::cout << Swarm_Method->_centers[0] << " " << Swarm_Method->_centers[1] << std::endl;
        EXPECT_TRUE(!_should_be_initialized);
        EXPECT_FALSE(!_should_not_be_initialized);
    }

    //Test for string update to work correctly
    
    //Create hypothetical drifts; add them to our initialized mock CVs, then check that each processor has the right values of reparameterized centers
    
    if(mpiid == 0)
        Swarm_Method->_cv_drift = {0.3, 1.0};
    else if(mpiid == 1)
        Swarm_Method->_cv_drift = {-0.4, 0.8};
    else if(mpiid == 2)
        Swarm_Method->_cv_drift = {-0.2, -0.5};
    
    //std::cout << " Hello 2 " << std::endl; //Debugging
    Swarm_Method->UpdateWorldString(cvlist_initialized);
    Swarm_Method->SetSendRecvNeighbors();
    Swarm_Method->StringUpdate();
    //std::cout << " Hello 3 " << std::endl; //Debugging

    //std::cout << mpiid << " " << Swarm_Method->_centers[0] << " " << Swarm_Method->_centers[1] << std::endl; //Debugging
    if(mpiid == 0)
    {
        EXPECT_NEAR(Swarm_Method->_centers[0], -0.19, eps);
        EXPECT_NEAR(Swarm_Method->_centers[1], 0.49, eps);
    }
    else if(mpiid == 1)
    {
        EXPECT_NEAR(Swarm_Method->_centers[0], 1.8882355, eps);
        EXPECT_NEAR(Swarm_Method->_centers[1], -0.966577, eps);
    }
    else if(mpiid == 2)
    {
        EXPECT_NEAR(Swarm_Method->_centers[0], 2.6, eps);
        EXPECT_NEAR(Swarm_Method->_centers[1], 1.7, eps);
    }

    std::cout.clear();
}

int main(int argc, char *argv[])
{
    int result = 0;
    ::testing::InitGoogleTest(&argc, argv);
    MPI_Init(&argc, &argv);
    result = RUN_ALL_TESTS();
    MPI_Finalize();


    return result;
}






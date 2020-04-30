#include "Tests.h"

#include "Methods/Swarm.h"
#include "CVs/MockCV.h"

using namespace SSAGES;

class SwarmTests : public ::testing::Test
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
        std::vector<double> centers_;
        centers_.resize(2); //Centers will contain two elements to check against
		if(mpiid == 0)
			centers_ = {-0.49, -0.51};
		else if(mpiid == 1)
			centers_ = {2.2, -1.8};
		else if(mpiid == 2)
			centers_ = {2.8, 2.2};
        std::vector<double> cvspring_;
        
        unsigned int iteration = 5;
		unsigned int numnodes = world.size();
		worldstring = {{-0.49, -0.51}, {2.2, -1.8}, {2.8, 2.2}};

        //Dummy Swarm constructor - world, comm, centers, maxiterations, cvspring, frequency, initial steps, harvest length, number trajectories, swarm length
        Swarm_Method = new Swarm(world, comm, centers_, 1000, cvspring_, 1, 10000, 10, 100, 20); 
        Swarm_Method->worldstring_ = worldstring;
        Swarm_Method->mpiid_ = mpiid;
        Swarm_Method->numnodes_ = numnodes;
        Swarm_Method->iteration_ = iteration;
    }

	virtual void TearDown() 
	{
		delete CV1;
		delete CV2;
		delete CV3;

        delete Swarm_Method;
	}

    mxx::comm world;
    mxx::comm comm = world.split(world.rank() < 3 ? world.rank() : 4);
    
    unsigned int mpiid;
    std::vector<std::vector<double>> worldstring;

	MockCV* CV1;
	MockCV* CV2;
	MockCV* CV3;

	CVList cvlist_initialized;
    CVList cvlist_not_initialized;
    
    Swarm* Swarm_Method;
};

TEST_F(SwarmTests,dummytest)
{
    //The intention here is to suppress output from functions
    std::cout.setstate(std::ios_base::failbit); 

	// only run this test with 3 ~faux walkers~
	ASSERT_TRUE(3 == world.size()) << "ERROR: Swarm unit test must be run with 3 processes";

    //std::cout << " Hello 1 " << std::endl; //Debugging
    //Test for CV Initialized
    if(mpiid == 0)
    {
        bool should_be_initialized = Swarm_Method->CVInitialized(cvlist_initialized);
        bool should_not_be_initialized = Swarm_Method->CVInitialized(cvlist_not_initialized);

        //Backwards logic could be cleaned up in swarm itself
        //std::cout << Swarm_Method->centers_[0] << " " << Swarm_Method->centers_[1] << std::endl;
        EXPECT_TRUE(should_be_initialized);
        EXPECT_FALSE(should_not_be_initialized);
    }

    //Test for string update to work correctly
    
    //Create hypothetical drifts; add them to our initialized mock CVs, then check that each processor has the right values of reparameterized centers
    
    if(mpiid == 0)
        Swarm_Method->cv_drift_ = {0.3, 1.0};
    else if(mpiid == 1)
        Swarm_Method->cv_drift_ = {-0.4, 0.8};
    else if(mpiid == 2)
        Swarm_Method->cv_drift_ = {-0.2, -0.5};
    
    //std::cout << " Hello 2 " << std::endl; //Debugging
    Swarm_Method->UpdateWorldString(cvlist_initialized);
    Swarm_Method->SetSendRecvNeighbors();
    Swarm_Method->StringUpdate();
    //std::cout << " Hello 3 " << std::endl; //Debugging

    //std::cout << mpiid << " " << Swarm_Method->centers_[0] << " " << Swarm_Method->centers_[1] << std::endl; //Debugging
    if(mpiid == 0)
    {
        EXPECT_NEAR(Swarm_Method->centers_[0], -0.19, eps);
        EXPECT_NEAR(Swarm_Method->centers_[1], 0.49, eps);
    }
    else if(mpiid == 1)
    {
        EXPECT_NEAR(Swarm_Method->centers_[0], 1.8882355, 1e4*eps);
        EXPECT_NEAR(Swarm_Method->centers_[1], -0.966577, 1e4*eps);
    }
    else if(mpiid == 2)
    {
        EXPECT_NEAR(Swarm_Method->centers_[0], 2.6, eps);
        EXPECT_NEAR(Swarm_Method->centers_[1], 1.7, eps);
    }

    std::cout.clear();
}

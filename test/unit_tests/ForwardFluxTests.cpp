#include "gtest/gtest.h"
#include <boost/mpi.hpp>

#include "../src/Snapshot.h"
#define private public
#define protected public
#include "../src/Methods/ForwardFlux.h"
#include "../src/CVs/MockCV.h"

using namespace SSAGES;

// Test calculation up to accuracy of 10^-10
const double eps = 0.0000000001;

class ForwardFluxTest : public ::testing::Test 
{
protected:

    virtual void SetUp() 
	{

		// Set up dummy CVs.
		//Vector3 grad1 = {1,1,0};

		//CV1 = new MockCV(-0.5,{1,0,0},-2,2); // In bounds.
		//CV2 = new MockCV(0.5,{0,1,0},-2,2); // In bounds.
		//CV3 = new MockCV(-2,{0,0,1},-2,2); // Out of bounds.

		//cvlist.push_back(CV1);
		//cvlist.push_back(CV2);
		//cvlist.push_back(CV3);	

		// Set up three timepoints.
		snapshot1 = new Snapshot(comm, 0);
		snapshot2 = new Snapshot(comm, 0);
		snapshot3 = new Snapshot(comm, 0);		

		snapshot1->SetNumAtoms(1);
		snapshot2->SetNumAtoms(1);
		snapshot3->SetNumAtoms(1);

		auto& loc1 = snapshot1 ->GetPositions();
		auto& loc2 = snapshot2 ->GetPositions();
		auto& loc3 = snapshot3 ->GetPositions();

		loc1.resize(1);
		loc2.resize(1);
		loc3.resize(1);

		loc1[0][0] = 2.0;
		loc1[0][1] = 3.0;
		loc1[0][2] = 4.0;
        //loc1[0] = {2,3,4}; //does this work?

		loc2[0][0] = 1.0;
		loc2[0][1] = 0.0;
		loc2[0][2] = -1.0;

		loc3[0][0] = 0.0;
		loc3[0][1] = -3.0;
		loc3[0][2] = -6.0;


		//forces1.resize(1);
		//forces1[0].resize(3);
		//forces1[0] = {0,0,0};

        //setup interfaces!!!
        // perhaps some that are increasing -1 -> +1
        // and also some that run +1 -> -1
        // currently this will break my FFS code
        fixme!



		// unsigned int mpiid1 = snapshot1 -> GetWalkerID();
		// unsigned int mpiid2 = snapshot2 -> GetWalkerID();
		// unsigned int mpiid3 = snapshot3 -> GetWalkerID();

		// mpiid1 = 0;
		// mpiid2 = 0;
		// mpiid3 = 0;

		// Set up ForwardFlux method.
		std::vector<std::vector<double>> histdetails;
		std::vector<std::vector<double>> restraint;
		std::vector<bool> isperiodic;
		std::vector<std::vector<double>> periodicboundaries;
		for(size_t i=0; i<cvlist.size(); ++i)
		{
			std::vector<double> temp1 = {-1, 1, 2};
			std::vector<double> temp2 = {-1.5, 1.5, 10};
			histdetails.push_back(temp1);
			restraint.push_back(temp2);
			isperiodic.push_back(false);
			periodicboundaries.push_back(std::vector<double>(2, 0.0));
		}

		Method = new ForwardFlux(world, // MPI global communicator
						 comm,  // MPI local communicator
						 restraint, // Min, max, spring constant
						 isperiodic, // Periodicity in CVs
						 periodicboundaries, // Position of periodic boundaries
						 5, // min hists before applying bias
						 false, // mass wheighing
						 "testout", // filename
						 histdetails, // hist details
						 10, // Backup interval
						 1, // Unit conversion
						 1, // timestep
						 1); // frequency
	}

	virtual void TearDown() 
	{
		delete snapshot1;
		delete snapshot2;
		delete snapshot3;

		//delete CV1;
		//delete CV2;
		//delete CV3;

		delete Method;
	}

	boost::mpi::communicator world;
	boost::mpi::communicator comm;

	Snapshot* snapshot1;
	Snapshot* snapshot2;
	Snapshot* snapshot3;

	MockCV* CV1;
	MockCV* CV2;
	MockCV* CV3;

    FFSConfigID ffsconfig1;
    FFSConfigID ffsconfig2;

	ForwardFlux* Method;	

	CVList cvlist;
};

TEST_F(ForwardFluxTest,FFSConfigID)
{
    //Check FFSConfigID

    //make sure constructor is working
    ffsconfig1 = FFSConfigID(1,2,3,4,5,6);
    EXPECT_NEAR(ffsconfig1.l,1, eps);
    EXPECT_NEAR(ffsconfig1.n,2, eps);
    EXPECT_NEAR(ffsconfig1.a,3, eps);
    EXPECT_NEAR(ffsconfig1.lprev,4, eps);
    EXPECT_NEAR(ffsconfig1.nprev,5, eps);
    EXPECT_NEAR(ffsconfig1.aprev,6, eps);
    
    // make sure operator equals is working
    ffsconfig2 = ffsconfig1;
    EXPECT_NEAR(ffsconfig2.l,1, eps);
    EXPECT_NEAR(ffsconfig2.n,2, eps);
    EXPECT_NEAR(ffsconfig2.a,3, eps);
    EXPECT_NEAR(ffsconfig2.lprev,4, eps);
    EXPECT_NEAR(ffsconfig2.nprev,5, eps);
    EXPECT_NEAR(ffsconfig2.aprev,6, eps);

}
TEST_F(ForwardFluxTest,CheckInitialStructure)
{
    //given a cv that is valid, make sure function agrees that its valid

    //given a cv that is invalid, make sure function agrees that its invalid
}

TEST_F(ForwardFluxTest,HasReturnedToA)
{
    //give A, expect true

    //give not A, expect false
}
TEST_F(ForwardFluxTest,HasCrossedInterface)
{
    //for interfaces increasing in the + direction
    //   set cvvalue and cvvalue_previous, expect 1
    //   set cvvalue and cvvalue_previous, expect 0
    //   set cvvalue and cvvalue_previous, expect -1

    //for interfaces increasing in the - direction
    //   set cvvalue and cvvalue_previous, expect 1
    //   set cvvalue and cvvalue_previous, expect 0
    //   set cvvalue and cvvalue_previous, expect -1
}

TEST_F(ForwardFluxTest,ReadWriteFFSConfiguration)
{
    //writeFFSconfig file
    //read it
    //check that values are the same

}

TEST_F(ForwardFluxTest,PopQueueMPI)
{
    //setup queue in MPI, SwarmTest is a good example for how to do this
    //pop it several times
    //make sure everyone has the correct configuration
}
TEST_F(ForwardFluxTest,ComputeCommittorProbability){
    //write some FFS files corresponding to some path

    //calculate the committor probabilities

    //check that they're the expected value
}

//--------------------------------------------------------------------
// SHOULD THESE TWO TESTS BE IN A DirectForwardFluxTests.cpp?
// For now I'll leave it here

TEST_F(ForwardFluxTest,InitializeQueue){
    //give Lambda0ConfigLibrary some values
    //seed  _generator
    // check that queue size
    // check the values of entries in queue

    //check that final PopQueue was called correctly by looking at myFFSConfigID

}

TEST_F(ForwardFluxTest,CheckForInterfaceCrossings){
    //call PostIntegration several times using ficticious snapshot files to 'evolve' the system forward in time
    //throughout this process, keep track of the queue, the number of successes, failures, current_interface, etc
    //this is probabily the most important unit test in the entire method

//--------------------------------------------------------------------


int main(int argc, char *argv[])
{
    ::testing::InitGoogleTest(&argc, argv);
     boost::mpi::environment env(argc,argv);
    int ret = RUN_ALL_TESTS();

    return ret;
}






#include "gtest/gtest.h"
#include <boost/mpi.hpp>

#define private public
#define protected public
#include "../src/Snapshot.h"
#include "../src/Methods/ForwardFlux.h"
#include "../src/Methods/DirectForwardFlux.h"
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
        // gradient isn't used in FFS so I set it to zero
		CV1 = new MockCV(-1,{0,0,0},-2,2); // In bounds.
		cvlist.push_back(CV1);

		// Set up three timepoints.
		snapshot1 = new Snapshot(comm, 0);
		snapshot1->SetNumAtoms(1);

		auto& loc1 = snapshot1 ->GetPositions();
		loc1.resize(1);
        loc1[0] = {2,3,4};
		auto& vel1 = snapshot1->GetVelocities();
        vel1.resize(1);
        vel1[0] = {5,6,7};
        snapshot1->_atomids.resize(1);
        snapshot1->_atomids[0] = 1;

		// unsigned int mpiid1 = snapshot1 -> GetWalkerID();
		// unsigned int mpiid2 = snapshot2 -> GetWalkerID();
		// unsigned int mpiid3 = snapshot3 -> GetWalkerID();

		// mpiid1 = 0;
		// mpiid2 = 0;
		// mpiid3 = 0;


		// Set up ForwardFlux method.
        unsigned int ninterfaces = 5;
        std::vector<double> interfaces;
        interfaces.resize(ninterfaces);
        interfaces[0]=-1.0;
        interfaces[1]=-0.95;
        interfaces[2]=-0.8;
        interfaces[3]= 0;
        interfaces[4]= 1.0;
        
        std::vector<unsigned int> M;
        M.resize(ninterfaces);
        for(unsigned int i=0;i<ninterfaces;i++) M[i] = 5;
        std::string output_directory="FFSoutput";

        mpirank = world.rank();

        //Method with increasing interfaces
		MethodA = new DirectForwardFlux(world, // MPI global communicator
						 comm,  // MPI local communicator
						 ninterfaces, // ninterfaces
						 interfaces, // interfaces vector
						 10, //N0target
						 M, // M vector
						 true, // initial flux flag
						 true, // savetrajectories
						 0, // currentinterface
						 output_directory, // output directory
						 1); // frequency

        //Method with decreasing interfaces
        interfaces[0]=1.0;
        interfaces[1]=0.95;
        interfaces[2]=0.8;
        interfaces[3]= 0;
        interfaces[4]=-1.0;
        MethodB = new DirectForwardFlux(world, // MPI global communicator
						 comm,  // MPI local communicator
						 ninterfaces, // ninterfaces
						 interfaces, // interfaces vector
						 10, //N0target
						 M, // M vector
						 true, // initial flux flag
						 true, // savetrajectories
						 0, // currentinterface
						 output_directory, // output directory
						 1); // frequency

	}

	virtual void TearDown() 
	{
		delete snapshot1;
		//delete snapshot2;
		//delete snapshot3;

		delete CV1;
		//delete CV3;

		delete MethodA;
		delete MethodB;

	}

	boost::mpi::communicator world;
	boost::mpi::communicator comm;

	Snapshot* snapshot1;
	//Snapshot* snapshot2;
	//Snapshot* snapshot3;

    unsigned int mpirank;

	MockCV* CV1;
	//MockCV* CV3;

    ForwardFlux::FFSConfigID ffsconfig1;
    ForwardFlux::FFSConfigID ffsconfig2;

	DirectForwardFlux* MethodA;	
	DirectForwardFlux* MethodB;	

	CVList cvlist;
};

TEST_F(ForwardFluxTest,FFSConfigID)
{
    //Check FFSConfigID

    //make sure constructor is working
    ffsconfig1 = ForwardFlux::FFSConfigID(1,2,3,4,5,6);
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
    //set the cv
    cvlist.clear();
    cvlist.push_back(new MockCV(-1.1,{0,0,0},-2,2)); //should succeed
    testing::internal::CaptureStdout();
    MethodA->CheckInitialStructure(cvlist);
    std::string outputValidCV = testing::internal::GetCapturedStdout();
    EXPECT_TRUE(outputValidCV == "Running initial Flux calculations\n");   

    //given a cv that is invalid, make sure function agrees that its invalid
    //set the cv
    cvlist.clear();
    cvlist.push_back(new MockCV(0,{0,0,0},-2,2)); //should fail and exit.
    EXPECT_EXIT(MethodA->CheckInitialStructure(cvlist), ::testing::ExitedWithCode(1), "Please provide an initial configuration in State A. Exiting ....\n");   
}

TEST_F(ForwardFluxTest,HasReturnedToA)
{
    //give A, expect true
    EXPECT_TRUE(MethodA->HasReturnedToA(-1.1));
    //give not A, expect false
    EXPECT_FALSE(MethodA->HasReturnedToA(-0.9));

    //give A, expect true
    EXPECT_TRUE(MethodB->HasReturnedToA(1.1));
    //give not A, expect false
    EXPECT_FALSE(MethodB->HasReturnedToA(0.9));

}

TEST_F(ForwardFluxTest,Constructor)
{
    // Check that various variables were correctly initialized in constructor
    EXPECT_TRUE(MethodA->_interfaces_increase);
    EXPECT_FALSE(MethodB->_interfaces_increase);

    EXPECT_TRUE(MethodA->_interfaces.size() == MethodA->_ninterfaces);
    EXPECT_TRUE(MethodA->_N.size() == MethodA->_ninterfaces);
    EXPECT_TRUE(MethodA->_A.size() == MethodA->_ninterfaces);
    EXPECT_TRUE(MethodA->_S.size() == MethodA->_ninterfaces);
    EXPECT_TRUE(MethodA->_P.size() == MethodA->_ninterfaces);
    
}



TEST_F(ForwardFluxTest,HasCrossedInterface)
{
    //for interfaces increasing in the + direction
    EXPECT_TRUE(MethodA->HasCrossedInterface(-0.94,-0.96,1) == 1);
    EXPECT_TRUE(MethodA->HasCrossedInterface(-0.96,-0.94,1) == -1);
    EXPECT_TRUE(MethodA->HasCrossedInterface(-0.94,-0.94,1) == 0);

    //for interfaces increasing in the - direction
    EXPECT_TRUE(MethodB->HasCrossedInterface(0.94,0.96,1) == 1);
    EXPECT_TRUE(MethodB->HasCrossedInterface(0.96,0.94,1) == -1);
    EXPECT_TRUE(MethodB->HasCrossedInterface(0.96,0.96,1) == 0);

}

TEST_F(ForwardFluxTest,ReadWriteFFSConfiguration)
{
    ffsconfig1 = ForwardFlux::FFSConfigID(1,2,3,4,5,6);
    //writeFFSconfig file
    MethodA->WriteFFSConfiguration(snapshot1,ffsconfig1,true);

    //change ffsconfig prev values (these should be overwritten in ReadFFSConfig)
    ffsconfig1.lprev = 77;
    ffsconfig1.nprev = 88;
    ffsconfig1.aprev = 99;

    //also change snapshot position (should be overwritten in ReadFFSConfig)
    auto& loc1 = snapshot1->GetPositions();
    loc1[0] = {555,666,777};
    auto& vel1 = snapshot1->GetVelocities();
    vel1[0] = {111,222,333};
   
    MPI_Barrier(world);

    //read it
    MethodA->ReadFFSConfiguration(snapshot1,ffsconfig1,true);

    //check that ReadFFSconfiguration overrode ffsconfig1.lprev
    EXPECT_TRUE(ffsconfig1.lprev == 4);
    EXPECT_TRUE(ffsconfig1.nprev == 5);
    EXPECT_TRUE(ffsconfig1.aprev == 6);
    
    //not make sure the positions and velocities were reloaded correclty
    EXPECT_NEAR(loc1[0][0],2,eps);
    EXPECT_NEAR(loc1[0][1],3,eps);
    EXPECT_NEAR(loc1[0][2],4,eps);
    EXPECT_NEAR(vel1[0][0],5,eps);
    EXPECT_NEAR(vel1[0][1],6,eps);
    EXPECT_NEAR(vel1[0][2],7,eps);

    //check that the failed trajectories can also be read and written correctly

}

TEST_F(ForwardFluxTest,PopQueueMPI)
{
    //setup queue in MPI
    
	//EXPECT_TRUE(3 == world.size()) << 
    if (world.size() != 3){
      std::cout<< "WARNING: ForwardFlux unit test should be run with 3 processes to test PopQueueMPI fully\n";
    }

    //write files, since PopQueue needs them
    ffsconfig1 = ForwardFlux::FFSConfigID(0,10,0,0,10,0);
    MethodA->WriteFFSConfiguration(snapshot1,ffsconfig1,true);
    ffsconfig1.n = 11;
    MethodA->WriteFFSConfiguration(snapshot1,ffsconfig1,true);
    ffsconfig1.n = 12;
    MethodA->WriteFFSConfiguration(snapshot1,ffsconfig1,true);
    ffsconfig1.n = 13;
    MethodA->WriteFFSConfiguration(snapshot1,ffsconfig1,true);

    //initialize queue
    MethodA->FFSConfigIDQueue.emplace_back(0,10,0,0,10,0);
    MethodA->FFSConfigIDQueue.emplace_back(0,11,0,0,1,0);
    MethodA->FFSConfigIDQueue.emplace_back(0,12,0,0,2,0);

    bool shouldpop_local = false; 
    //pop rank 1 and 3
    if (mpirank == 0) shouldpop_local = true;
    else if (mpirank == 1){
      MethodA->myFFSConfigID = ForwardFlux::FFSConfigID(0,1,0,0,1,0);
      shouldpop_local = false;
    }
    else if (mpirank == 2) shouldpop_local = true;

    MethodA->PopQueueMPI(snapshot1, cvlist, shouldpop_local);
    
    //make sure everyone has the correct ffsID
    if (mpirank == 0)
        EXPECT_TRUE(MethodA->myFFSConfigID.n == 10);
    else if (mpirank == 1)
        EXPECT_TRUE(MethodA->myFFSConfigID.n == 1);
    else if (mpirank == 2)
        EXPECT_TRUE(MethodA->myFFSConfigID.n == 11);

    //make sure everyone has correct queue, depending on world size
    if (world.size() == 3){
      EXPECT_TRUE(MethodA->FFSConfigIDQueue.front().n == 12);
      EXPECT_TRUE(MethodA->FFSConfigIDQueue.back().n == 12);
    }
    else if (world.size() == 1){
      EXPECT_TRUE(MethodA->FFSConfigIDQueue.front().n == 11);
      EXPECT_TRUE(MethodA->FFSConfigIDQueue.back().n == 12);
    }
    
    //make sure trajectory file is open, if queue was popped
    if (shouldpop_local)
      EXPECT_TRUE(MethodA->_trajectory_file.is_open());

    //make sure cvprev == cv
    EXPECT_TRUE(MethodA->_cvvalue_previous == MethodA->_cvvalue);

    //now empty queue completely, make sure that _pop_tried_but_empty_queue is set correclty
    if (world.size() == 3){
      MethodA->_trajectory_file.close();
      MethodA->PopQueueMPI(snapshot1, cvlist, shouldpop_local);
      if ((mpirank == 0) || (mpirank == 1)){
        EXPECT_FALSE(MethodA->_pop_tried_but_empty_queue);
      }
      else if (mpirank == 2){
        EXPECT_TRUE(MethodA->_pop_tried_but_empty_queue);
      }
    }
    else if (world.size() == 1){
      MethodA->_trajectory_file.close();
      MethodA->PopQueueMPI(snapshot1, cvlist, shouldpop_local);
      MethodA->_trajectory_file.close();
      MethodA->PopQueueMPI(snapshot1, cvlist, shouldpop_local);
      EXPECT_FALSE(MethodA->_pop_tried_but_empty_queue);
      MethodA->PopQueueMPI(snapshot1, cvlist, shouldpop_local);
      EXPECT_TRUE(MethodA->_pop_tried_but_empty_queue);
      MethodA->_trajectory_file.close();
    }


}
TEST_F(ForwardFluxTest,ComputeCommittorProbability)
{
    //write some FFS files corresponding to some path

    //calculate the committor probabilities

    //check that they're the expected value
}

//--------------------------------------------------------------------
// SHOULD THESE TWO TESTS BE IN A DirectForwardFluxTests.cpp?
// For now I'll leave it here

TEST_F(ForwardFluxTest,InitializeQueue)
{
    //give Lambda0ConfigLibrary some values
    //seed  _generator
    // check that queue size
    // check the values of entries in queue

    //check that final PopQueue was called correctly by looking at myFFSConfigID

}

TEST_F(ForwardFluxTest,CheckForInterfaceCrossings)
{
    //call PostIntegration several times using ficticious snapshot files to 'evolve' the system forward in time
    //throughout this process, keep track of the queue, the number of successes, failures, current_interface, etc
    //this is probabily the most important unit test in the entire method

    //here I perform 3 "steps" or calls to CheckForInterfaceCrossings
    // each time, different logic should be triggered

  if (world.size() == 3){
    unsigned current_interface = 2;
    MethodA->_current_interface = current_interface;
    MethodA->_A[current_interface] = 0;
    MethodA->_S[current_interface] = 0;
    //M[2] = 2;

    ffsconfig1 = ForwardFlux::FFSConfigID(2,3,0,1,0,0);
    MethodA->WriteFFSConfiguration(snapshot1,ffsconfig1,true);
    ffsconfig1 = ForwardFlux::FFSConfigID(2,4,0,1,6,0);
    MethodA->WriteFFSConfiguration(snapshot1,ffsconfig1,true);
    //ffsconfig1 = ForwardFlux::FFSConfigID(2,0,1,1,7,0);
    //MethodA->WriteFFSConfiguration(snapshot1,ffsconfig1,true);

    //write files, these are the two successes
    //ffsconfig1 = ForwardFlux::FFSConfigID(3,0,0,2,9,0);
    //MethodA->WriteFFSConfiguration(snapshot1,ffsconfig1,true);
    //ffsconfig1 = ForwardFlux::FFSConfigID(3,1,0,2,10,0);
    //MethodA->WriteFFSConfiguration(snapshot1,ffsconfig1,true);

    if (mpirank == 0) MethodA->myFFSConfigID = ForwardFlux::FFSConfigID(2,0,0,1,2,0);
    if (mpirank == 1) MethodA->myFFSConfigID = ForwardFlux::FFSConfigID(2,1,0,1,3,0);
    if (mpirank == 2) MethodA->myFFSConfigID = ForwardFlux::FFSConfigID(2,2,0,1,4,0);

    //initialize queue
    //once the queue is empty, this interface is done
    MethodA->FFSConfigIDQueue.emplace_back(2,3,0,1,5,0);
    MethodA->FFSConfigIDQueue.emplace_back(2,4,0,1,6,0); 
    //MethodA->FFSConfigIDQueue.emplace_back(2,0,1,1,7,0);
   

   //------------------------------
   // Step 1
   //------------------------------
    //set cvprevious 
    if (mpirank == 0) MethodA->_cvvalue_previous = -0.1; //will succeed
    if (mpirank == 1) MethodA->_cvvalue_previous = -0.9; //will fail
    if (mpirank == 2) MethodA->_cvvalue_previous = -0.5; //will do neither
    
    //set cv's then run a step
    cvlist.clear();
    //interface 3 is at zero
    if (mpirank == 0) cvlist.push_back(new MockCV(0.1,{0,0,0},-2,2)); //should succeed
    if (mpirank == 1) cvlist.push_back(new MockCV(-0.9,{0,0,0},-2,2)); //should do neither
    if (mpirank == 2) cvlist.push_back(new MockCV(-0.5,{0,0,0},-2,2)); //snould do neither

    MethodA->CheckForInterfaceCrossings(snapshot1,cvlist);
    
    //check everyones successes and attempts count
    EXPECT_NEAR(MethodA->_current_interface,2,eps);    //interface is same
    EXPECT_NEAR(MethodA->_S[current_interface],1,eps); //successes inc by one
    EXPECT_NEAR(MethodA->_A[current_interface],1,eps); //attempts inc by one

    //check that each proc has the corrects FFSConfigID
    if (mpirank == 0){
      EXPECT_NEAR(MethodA->myFFSConfigID.l, 2, eps); //this is new off queue
      EXPECT_NEAR(MethodA->myFFSConfigID.n, 3, eps);
      EXPECT_NEAR(MethodA->myFFSConfigID.a, 0, eps);
    }
    else if (mpirank == 1){
      EXPECT_NEAR(MethodA->myFFSConfigID.l, 2, eps); //still previous ffsconfig
      EXPECT_NEAR(MethodA->myFFSConfigID.n, 1, eps);
      EXPECT_NEAR(MethodA->myFFSConfigID.a, 0, eps);
    }
    else if (mpirank == 2){
      EXPECT_NEAR(MethodA->myFFSConfigID.l, 2, eps); //still previous ffsconfig
      EXPECT_NEAR(MethodA->myFFSConfigID.n, 2, eps);
      EXPECT_NEAR(MethodA->myFFSConfigID.a, 0, eps);
    }


   //------------------------------
   // Step 2
   //------------------------------

    //this is an important next step, need to have 2 more queue pops
    // the first one will get the last config in the queue, the 2nd will have the _popqueue empty set

    //again set cvprevious (I'll use same as above for ease)
    if (mpirank == 0) MethodA->_cvvalue_previous = -0.1; //will succeed
    if (mpirank == 1) MethodA->_cvvalue_previous = -0.9; //will fail
    if (mpirank == 2) MethodA->_cvvalue_previous = -0.5; //will do neither
    //now set current cv
    cvlist.clear();
    if (mpirank == 0) cvlist.push_back(new MockCV(0.1,{0,0,0},-2,2)); //should succeed
    if (mpirank == 1) cvlist.push_back(new MockCV(-1.1,{0,0,0},-2,2)); //should fail
    if (mpirank == 2) cvlist.push_back(new MockCV(-0.5,{0,0,0},-2,2)); //snould do neither

    MethodA->CheckForInterfaceCrossings(snapshot1,cvlist);
    

    //check everyones successes and attempts count
    EXPECT_NEAR(MethodA->_current_interface,2,eps);    //interface is same
    EXPECT_NEAR(MethodA->_S[current_interface],2,eps); //successes inc by one
    EXPECT_NEAR(MethodA->_A[current_interface],3,eps); //attempts inc by two

    if (mpirank == 0){
      EXPECT_NEAR(MethodA->myFFSConfigID.l, 2, eps); //this is new off queue
      EXPECT_NEAR(MethodA->myFFSConfigID.n, 4, eps);
      EXPECT_NEAR(MethodA->myFFSConfigID.a, 0, eps);
    }
    else if (mpirank == 1){
      EXPECT_NEAR(MethodA->myFFSConfigID.l, 2, eps); //this is a prev config, but is a zombie job
      EXPECT_NEAR(MethodA->myFFSConfigID.n, 1, eps);
      EXPECT_NEAR(MethodA->myFFSConfigID.a, 0, eps);
      EXPECT_TRUE(MethodA->_pop_tried_but_empty_queue);//this is check for zombie job
    }
    else if (mpirank == 2){
      EXPECT_NEAR(MethodA->myFFSConfigID.l, 2, eps); //this is still the original, since so success or failure
      EXPECT_NEAR(MethodA->myFFSConfigID.n, 2, eps);
      EXPECT_NEAR(MethodA->myFFSConfigID.a, 0, eps);
    }
    //queue should be empty
    EXPECT_TRUE(MethodA->FFSConfigIDQueue.size()==0);

   //------------------------------
   // Step 3
   //------------------------------
    //then  have the last two jobs finish (on proc 0 and 2), this interfaces should be done

    //again set cvprevious 
    if (mpirank == 0) MethodA->_cvvalue_previous = -0.9; //will fail
    if (mpirank == 1) MethodA->_cvvalue_previous = -0.1; //is zombie, but I'll have it succeed, just to make sure its ignored
    if (mpirank == 2) MethodA->_cvvalue_previous = -0.1; //will succeed
    //now set current cv
    cvlist.clear();
    if (mpirank == 0) cvlist.push_back(new MockCV(-1.1,{0,0,0},-2,2)); //should fail
    if (mpirank == 1) cvlist.push_back(new MockCV( 0.1,{0,0,0},-2,2)); //should do nothing
    if (mpirank == 2) cvlist.push_back(new MockCV( 0.1,{0,0,0},-2,2)); //snould succeed

    MethodA->CheckForInterfaceCrossings(snapshot1,cvlist);
    current_interface = MethodA->_current_interface; //need to store in local var
 
    //check everyones successes and attempts count
    EXPECT_NEAR(MethodA->_current_interface,3,eps);    //interface is now +1
    EXPECT_NEAR(MethodA->_S[current_interface-1],3,eps); //successes of prev interface inc by one
    EXPECT_NEAR(MethodA->_A[current_interface-1],5,eps); //attempts of prev interface inc by two
    EXPECT_NEAR(MethodA->_S[current_interface],0,eps); //successes of current interfaces 0
    EXPECT_NEAR(MethodA->_A[current_interface],0,eps); //attempts of current interface 0

    //queue should have 5 new jobs, but 3 were popped so should have size 2
   EXPECT_NEAR(MethodA->FFSConfigIDQueue.size(),2,eps);

    // want to check on individual jobs?, tough because config were put in queue randomly
    if (mpirank == 0){

    }
    else if (mpirank == 1){
      EXPECT_FALSE(MethodA->_pop_tried_but_empty_queue);//make sure no longer zombie
    }
    else if (mpirank == 2){

    }


 }

}

//--------------------------------------------------------------------



int main(int argc, char *argv[])
{
    ::testing::InitGoogleTest(&argc, argv);
     boost::mpi::environment env(argc,argv);
    int ret = RUN_ALL_TESTS();

    return ret;
}






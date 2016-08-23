#include "gtest/gtest.h"
#include <boost/mpi.hpp>

#include "../src/Snapshot.h"


#define private public
#include "../src/Methods/ABF.h"
#include "../src/CVs/MockCV.h"

using namespace SSAGES;

// Test calculation up to accuracy of 10^-10
const double eps = 0.0000000001;

class ABFTest : public ::testing::Test 
{

protected:
    	virtual void SetUp() 
	{

	// Set up dummy CVs.
	Vector3 grad1 = {1,1,0};

	CV1 = new MockCV(-0.5,{1,0,0},-2,2); // In bounds.
	CV2 = new MockCV(0.5,{0,1,0},-2,2); // In bounds.
	CV3 = new MockCV(-2,{0,0,1},-2,2); // Out of bounds.
	
	cvlist.push_back(CV1);
	cvlist.push_back(CV2);
	cvlist.push_back(CV3);	

	// Set up three timepoints.
	snapshot1 = new Snapshot(comm, 0);
	snapshot2 = new Snapshot(comm, 0);
	snapshot3 = new Snapshot(comm, 0);		

	auto& loc1 = snapshot1 ->GetPositions();
	auto& loc2 = snapshot1 ->GetPositions();
	auto& loc3 = snapshot1 ->GetPositions();

	loc1.resize(1);
	loc2.resize(1);
	loc3.resize(1);

	auto& vel1 = snapshot1 -> GetVelocities();
	auto& vel2 = snapshot2 -> GetVelocities();
	auto& vel3 = snapshot3 -> GetVelocities();

	vel1.resize(1);
	vel2.resize(1);
	vel3.resize(1);
	
	vel1[0][0] = 2.0;
	vel1[0][1] = 3.0;
	vel1[0][2] = 4.0;

	vel2[0][0] = 1.0;
	vel2[0][1] = 0.0;
	vel2[0][2] = -1.0;

	vel3[0][0] = 0.0;
	vel3[0][1] = -3.0;
	vel3[0][2] = -6.0;

	auto& forces1 = snapshot1 -> GetForces();
	auto& forces2 = snapshot2 -> GetForces();
	auto& forces3 = snapshot3 -> GetForces();

	forces1.resize(1);
		forces1[0].resize(3);
	forces2.resize(1);
		forces2[0].resize(3);
	forces3.resize(1);
		forces3[0].resize(3);

	auto& mass1 = snapshot1 -> GetMasses();
	auto& mass2 = snapshot2 -> GetMasses();
	auto& mass3 = snapshot3 -> GetMasses();

	mass1.resize(1);
	mass2.resize(1);
	mass3.resize(1);

	mass1[0] = 2.0;
	mass2[0] = 2.0;
	mass3[0] = 2.0;

	// unsigned int mpiid1 = snapshot1 -> GetWalkerID();
	// unsigned int mpiid2 = snapshot2 -> GetWalkerID();
	// unsigned int mpiid3 = snapshot3 -> GetWalkerID();

	// mpiid1 = 0;
	// mpiid2 = 0;
	// mpiid3 = 0;


	// Set up ABF method.
	std::vector<std::vector<double>> histdetails;
	std::vector<std::vector<double>> restraint;
	std::vector<int> printdetails = {10, 0, 0, 0, 0, 0, 0, 0, 0};
	for(size_t i=0; i<cvlist.size(); ++i)
		{
		std::vector<double> temp1 = {-1, 1, 2};
		std::vector<double> temp2 = {-1.5, 1.5, 10};
		histdetails.push_back(temp1);
		restraint.push_back(temp2);
		}
	
		
	Method = new ABF(world, comm, histdetails, restraint, 1, 5, "testout", printdetails, 10, 1, 0, 1);	

	}

	virtual void TearDown() 
	{
	delete snapshot1;
	delete snapshot2;
	delete snapshot3;

	delete CV1;
	delete CV2;
	delete CV3;

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

	ABF* Method;	

	CVList cvlist;

};

TEST_F(ABFTest,dummytest)
{
	// Initialize CVs.	
	CV1->MockCV::Initialize(*snapshot1);
	CV2->MockCV::Initialize(*snapshot1);
	CV3->MockCV::Initialize(*snapshot1);

	Method->ABF::PreSimulation(snapshot1, cvlist);

	// Test dimensionality is correct.
	EXPECT_TRUE(Method->_dim == 3);
	int dim = Method->_dim;

	// Check for correct initialization.
	EXPECT_TRUE(Method->_N.size() == 8);
	EXPECT_TRUE(Method->_F.size() == 24);

	for(size_t i = 0; i < (Method->_N).size() ; ++i)
		{
		EXPECT_TRUE(Method->_N[i] == 0);

		EXPECT_TRUE(Method->_F[dim*i] == 0);
		EXPECT_TRUE(Method->_F[dim*i+1] == 0);
		EXPECT_TRUE(Method->_F[dim*i+2] == 0);
		}

	// Check Setters
	Method->SetHistogram({1,1},{1,1});

	EXPECT_TRUE(Method->_N[0] == 1);
	EXPECT_TRUE(Method->_N[1] == 1);
	EXPECT_TRUE(Method->_F[0] == 1);
	EXPECT_TRUE(Method->_F[1] == 1);

	// Reset
	Method->_N.resize(0);
	Method->_F.resize(0);
	Method->ABF::PreSimulation(snapshot1, cvlist);

	// Check for correct re-initialization.
	EXPECT_TRUE(Method->_N.size() == 8);
	EXPECT_TRUE(Method->_F.size() == 24);

	for(size_t i = 0; i < (Method->_N).size() ; ++i)
		{
		EXPECT_TRUE(Method->_N[i] == 0);

		EXPECT_TRUE(Method->_F[dim*i] == 0);
		EXPECT_TRUE(Method->_F[dim*i+1] == 0);
		EXPECT_TRUE(Method->_F[dim*i+2] == 0);
		}
	
	// Take one MD step.
	Method->ABF::PostIntegration(snapshot1, cvlist);

	// Test whether bound checking works. There should be no hits as 3rd CV was out of bounds.
	for(size_t i = 0; i < (Method->_N).size() ; ++i)
		{
		EXPECT_TRUE(Method->_N[i] == 0);

		EXPECT_TRUE(Method->_F[dim*i] == 0);
		EXPECT_TRUE(Method->_F[dim*i+1] == 0);
		EXPECT_TRUE(Method->_F[dim*i+2] == 0);
		}

	// The only bias force should be from CV restraints.
	double bias = (CV3->GetValue()+1.5)*10;

	for(size_t i = 0; i < Method->_biases[0].size() ; ++i)
		{
		EXPECT_TRUE(-Method->_biases[0][i] == bias*CV3->GetGradient()[0][i]);
		}

	// Put CVs in bound.
	CV3->_val = -1;
		
	// Take two more MD steps
	Method->ABF::PostIntegration(snapshot2, cvlist);
	Method->ABF::PostIntegration(snapshot3, cvlist);

	// Check that hits are recorded correctly
	EXPECT_TRUE(Method->_N[2] == 2);

	EXPECT_TRUE(Method->_F[2*dim] == -8.0);
	EXPECT_TRUE(Method->_F[2*dim+1] == -20.4);
	EXPECT_TRUE(Method->_F[2*dim+2] == -32.8);

	// Check that the biases are accurate
	EXPECT_TRUE(Method->_biases[0][0] == 1.6);
	EXPECT_TRUE(Method->_biases[0][1] == 4.08);
	EXPECT_TRUE(Method->_biases[0][2] == 6.56);

	// Check that the biases added to forces correctly
	EXPECT_TRUE(snapshot3->GetForces()[0][0] == 1.6);
	EXPECT_TRUE(snapshot3->GetForces()[0][1] == 4.08);
	EXPECT_TRUE(snapshot3->GetForces()[0][2] == 6.56);

}

int main(int argc, char *argv[])
{
    ::testing::InitGoogleTest(&argc, argv);
     boost::mpi::environment env(argc,argv);
    int ret = RUN_ALL_TESTS();

    return ret;
}






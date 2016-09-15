#include "gtest/gtest.h"
#include <boost/mpi.hpp>

#include "../src/Snapshot.h"


#define private public
#define protected public 

#include "../src/Methods/FiniteTempString.h"
#include "../src/CVs/MockCV.h"

using namespace SSAGES;

/* Test reparametrization up to accuracy of 0.3:
 * Larger than normal eps since there is some discrepancy between
 * different cubic spline interpolators */
const double reparam_eps = 0.3;

class FTSTest : public ::testing::Test 
{
protected:

    virtual void SetUp() 
	{

		// Set up dummy CVs.
		CV1 = new MockCV(-0.5,{1,0,0},-2,2); 
		CV2 = new MockCV(-0.5,{0,1,0},-2,2);

		cvlist.push_back(CV1);
		cvlist.push_back(CV2);

        mpiid = world.rank();

        // Set up FTS parameters.	
        _centers.resize(2); // Centers will contain two elements to check against
		if(mpiid == 0)
			_centers = {-0.49, -0.51};
		else if(mpiid == 1)
			_centers = {2.2, -1.8};
		else if(mpiid == 2)
			_centers = {2.8, 2.2};
		else if(mpiid == 3)
			_centers = {4.5, 3.7};

		worldstring = {{-0.49, -0.51}, {2.2, -1.8}, {2.8, 2.2}, {4.5, 3.7}};

		/* FTS dummy constructor: 	world, comm, centers, maxiterations, blockiterations, tau, 
		 * 							cvspring, kappa, spring iter, frequency */
		FTS_Method = new FiniteTempString(world, comm, _centers, 5, 100, 0.1, _cvspring, 0.1, 1, 1);
		FTS_Method->_worldstring = worldstring;
		FTS_Method->_mpiid = mpiid;
		FTS_Method->_numnodes = world.size();
    }

	virtual void TearDown() 
	{
		delete CV1;
		delete CV2;

		delete FTS_Method;
	}

    boost::mpi::communicator world;
    boost::mpi::communicator comm = world.split(world.rank() < 4 ? world.rank() : 5);
    
    unsigned int mpiid;
	std::vector<double> _centers;
	std::vector<double> _cvspring;
	std::vector<std::vector<double>> worldstring;

	MockCV* CV1;
	MockCV* CV2;

	CVList cvlist;
    
	FiniteTempString* FTS_Method;
};

TEST_F(FTSTest,dummytest)
{
    //The intention here is to suppress output from functions
    std::cout.setstate(std::ios_base::failbit); 

	/* Only run this test with 4 ~faux walkers~
	 * Other string method unit tests use 3; using 4 here since reparametrization
	 * values are being compared with those generated using scipy.interpolate.interp1d,
	 * which requires at least 4 data points to perform a cubic spline */
	ASSERT_TRUE(4 == world.size()) << "ERROR: FTS unit test must be run with 4 processes";

	/* Test FTS_Method->InCell(cvlist)
	 * Check if current CVs is within current Voronoi cell, as defined by centers in _worldstring */
	if(mpiid == 0){
		EXPECT_TRUE(FTS_Method->InCell(cvlist)) << "ERROR: CVs are within cell but InCell() returns false";
	} else {
		EXPECT_FALSE(FTS_Method->InCell(cvlist)) << "ERROR: CVs are out of cell but InCell() returns true";
	}

	/* Test FTS_Method->StringUpdate()
	 * Test that string is updated correctly after moving toward running averages and 
	 * equal arc length reparametrization */
	if(mpiid == 0)
		FTS_Method->_newcenters = {-1, -1};
	else if (mpiid == 1)
		FTS_Method->_newcenters = {0, 0};
	else if (mpiid == 2)
		FTS_Method->_newcenters = {3, 2};
	else if (mpiid == 3)
		FTS_Method->_newcenters = {5, 4};

	FTS_Method->SetSendRecvNeighbors();
	FTS_Method->StringUpdate();

	if(mpiid == 0){
		EXPECT_NEAR(FTS_Method->_centers[0], -0.541, reparam_eps) << "ERROR: StringUpdate() updates _centers incorrectly";
		EXPECT_NEAR(FTS_Method->_centers[1], -0.559, reparam_eps) << "ERROR: StringUpdate() updates _centers incorrectly";
	} else if (mpiid == 1){
		EXPECT_NEAR(FTS_Method->_centers[0], 1.9631, reparam_eps) << "ERROR: StringUpdate() updates _centers incorrectly";
		EXPECT_NEAR(FTS_Method->_centers[1], -1.3254, reparam_eps) << "ERROR: StringUpdate() updates _centers incorrectly";
	} else if (mpiid == 2){
		EXPECT_NEAR(FTS_Method->_centers[0], 2.6365, reparam_eps) << "ERROR: StringUpdate() updates _centers incorrectly";
		EXPECT_NEAR(FTS_Method->_centers[1], 1.2956, reparam_eps) << "ERROR: StringUpdate() updates _centers incorrectly";
	} else if (mpiid == 3){
		EXPECT_NEAR(FTS_Method->_centers[0], 4.2824, reparam_eps) << "ERROR: StringUpdate() updates _centers incorrectly";
		EXPECT_NEAR(FTS_Method->_centers[1], 3.5016, reparam_eps) << "ERROR: StringUpdate() updates _centers incorrectly";
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






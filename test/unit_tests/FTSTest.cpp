#include "Tests.h"

#include "CVs/CVManager.h"
#include "Methods/FiniteTempString.h"
#include "CVs/MockCV.h"

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

		cvmanager.AddCV(CV1);
		cvmanager.AddCV(CV2);

        mpiid = world.rank();

        // Set up FTS parameters.	
        centers_.resize(2); // Centers will contain two elements to check against
		if(mpiid == 0)
			centers_ = {-0.49, -0.51};
		else if(mpiid == 1)
			centers_ = {2.2, -1.8};
		else if(mpiid == 2)
			centers_ = {2.8, 2.2};
		else if(mpiid == 3)
			centers_ = {4.5, 3.7};

		worldstring = {{-0.49, -0.51}, {2.2, -1.8}, {2.8, 2.2}, {4.5, 3.7}};

		/* FTS dummy constructor: 	world, comm, centers, maxiterations, blockiterations, tau, 
		 * 							cvspring, kappa, spring iter, frequency */
		FTS_Method = new FiniteTempString(world, comm, centers_, 5, 100, 0.1, cvspring_, 0.1, 1, 1);
		FTS_Method->worldstring_ = worldstring;
		FTS_Method->mpiid_ = mpiid;
		FTS_Method->numnodes_ = world.size();
    }

	virtual void TearDown() 
	{
		delete FTS_Method;
	}

    mxx::comm world;
    mxx::comm comm = world.split(world.rank() < 4 ? world.rank() : 5);
    
    unsigned int mpiid;
	std::vector<double> centers_;
	std::vector<double> cvspring_;
	std::vector<std::vector<double>> worldstring;

	MockCV* CV1;
	MockCV* CV2;

	CVManager cvmanager;
    
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
	 * Check if current CVs is within current Voronoi cell, as defined by centers in worldstring_ */
	if(mpiid == 0){
		EXPECT_TRUE(FTS_Method->InCell(cvmanager.GetCVs())) << "ERROR: CVs are within cell but InCell() returns false";
	} else {
		EXPECT_FALSE(FTS_Method->InCell(cvmanager.GetCVs())) << "ERROR: CVs are out of cell but InCell() returns true";
	}

	/* Test FTS_Method->StringUpdate()
	 * Test that string is updated correctly after moving toward running averages and 
	 * equal arc length reparametrization */
	if(mpiid == 0)
		FTS_Method->newcenters_ = {-1, -1};
	else if (mpiid == 1)
		FTS_Method->newcenters_ = {0, 0};
	else if (mpiid == 2)
		FTS_Method->newcenters_ = {3, 2};
	else if (mpiid == 3)
		FTS_Method->newcenters_ = {5, 4};

	FTS_Method->SetSendRecvNeighbors();
	FTS_Method->StringUpdate();

	if(mpiid == 0){
		EXPECT_NEAR(FTS_Method->centers_[0], -0.541, reparam_eps) << "ERROR: StringUpdate() updates centers_ incorrectly";
		EXPECT_NEAR(FTS_Method->centers_[1], -0.559, reparam_eps) << "ERROR: StringUpdate() updates centers_ incorrectly";
	} else if (mpiid == 1){
		EXPECT_NEAR(FTS_Method->centers_[0], 1.9631, reparam_eps) << "ERROR: StringUpdate() updates centers_ incorrectly";
		EXPECT_NEAR(FTS_Method->centers_[1], -1.3254, reparam_eps) << "ERROR: StringUpdate() updates centers_ incorrectly";
	} else if (mpiid == 2){
		EXPECT_NEAR(FTS_Method->centers_[0], 2.6365, reparam_eps) << "ERROR: StringUpdate() updates centers_ incorrectly";
		EXPECT_NEAR(FTS_Method->centers_[1], 1.2956, reparam_eps) << "ERROR: StringUpdate() updates centers_ incorrectly";
	} else if (mpiid == 3){
		EXPECT_NEAR(FTS_Method->centers_[0], 4.2824, reparam_eps) << "ERROR: StringUpdate() updates centers_ incorrectly";
		EXPECT_NEAR(FTS_Method->centers_[1], 3.5016, reparam_eps) << "ERROR: StringUpdate() updates centers_ incorrectly";
	}
    
	std::cout.clear();
}

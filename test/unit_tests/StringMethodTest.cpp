#include "gtest/gtest.h"
#include <boost/mpi.hpp>

#define private public
#define protected public 

#include "../src/Methods/FiniteTempString.h"
#include "../src/Methods/Swarm.h"
#include "../src/CVs/MockCV.h"

using namespace SSAGES;

// Test calculation up to accuracy of 10^-5
const double eps = 0.00001;
const double sq_dist_eps = 0.001; //Is this too big?  We might not actually get accurate enough results out of square dist
class StringMethodTest : public ::testing::Test
{
	protected:

		virtual void SetUp(){
			// Vectors for testing distance function
			vec1 = {-1, 0, 1};
			vec2 = {-2, 3, 4};

			// Tolerance and max iterations 
			maxiterator = 5;
			iteration = 5;

			// Set up dummy snapshot
			snapshot1 = new Snapshot(comm, 0);

			// Centers & worldstring
			numnodes = world.size();
	    	mpiid = world.rank();
			worldstring = {{-1, -1}, {2, -2}, {3, 2}};

			if(mpiid == 0)
				_centers = {-1.1, -1.1};
			else if(mpiid == 1)
				_centers = {2.2, -1.8};
			else if(mpiid == 2)
				_centers = {2.8, 2.2};

			// Set up mock CVs
			CV1 = new MockCV(-0.5, {1, 0, 0}, -2, 2); 
			CV2 = new MockCV(0.5, {0, 1, 0}, -2, 2);
			cvlist.push_back(CV1);
			cvlist.push_back(CV2);
        	
			/* Need a dummy method to test StringMethods functions because we can't create an instance 
			 * of a class with purely virtual functions (as StringMethod has).  
			 * Prefer FTS for consistency when calling common functions. */

        	// Dummy FTS constructor - world, comm, centers, maxiterations, block iterations, tau, cv spring, kappa, spring iter, frequency 
        	FTS_Method = new FiniteTempString(world, comm, _centers, maxiterator, 100, 1, _cvspring, 1, 1, 1);
        	FTS_Method->_worldstring = worldstring;
        	FTS_Method->_mpiid = mpiid;
        	FTS_Method->_numnodes = numnodes;
        	FTS_Method->_iteration = iteration;
		}

		virtual void TearDown(){
			delete CV1;
			delete CV2;
            
            delete snapshot1;

            delete FTS_Method;
		}

		boost::mpi::communicator world;
		boost::mpi::communicator comm = world.split(world.rank() < 3 ? world.rank() : 4);
        
        std::vector<double> vec1, vec2, _centers, _cvspring;
		std::vector<std::vector<double>> worldstring;
		unsigned int mpiid, maxiterator, iteration, numnodes;

    	MockCV* CV1;
	    MockCV* CV2;

		CVList cvlist;
        
        Snapshot* snapshot1;
	
        FiniteTempString* FTS_Method;
};

TEST_F(StringMethodTest,dummytest)
{
    //The intention here is to suppress output from functions
    std::cout.setstate(std::ios_base::failbit); 

	// only run this test with 3 ~faux walkers~
	ASSERT_TRUE(3 == world.size()) << "ERROR: String Method unit test must be run with 3 processes" ;

	// Test StringMethod->sqdist()
	ASSERT_EQ(vec1.size(), vec2.size());
	EXPECT_NEAR(FTS_Method->sqdist(vec1, vec2), 4.359, sq_dist_eps);
	// currently stringmethod::sqdist doesnt check if its two inputs are different sizes
	//
	// also: somewhere along the way this function stopped returning the squared distance
	// and returns [the already square-rooted] distance... function name is confusing
	//
	// eigen???
	
	// Test StringMethod->TolCheck()
	// Default is no tolerance defined: _tol is empty
	EXPECT_FALSE(FTS_Method->TolCheck()) << "ERROR: TolCheck() fails if no _tol value provided";
	// Initialize tolerance values: begin with tolerance criteria not met 
    FTS_Method->_tol = {0.0001, 0.0001};
	EXPECT_FALSE(FTS_Method->TolCheck()) << "ERROR: TolCheck() gives false positive";
	// Only one CV within tolerance:
	FTS_Method->_tol[0] = 10;
	EXPECT_FALSE(FTS_Method->TolCheck()) << "ERROR: TolCheck() gives false positive";
	FTS_Method->_tol[1] = 10;
	// After tolerance criteria is met:
	EXPECT_TRUE(FTS_Method->TolCheck()) << "ERROR: TolCheck() gives false negative";
	
	// Test StringMethod->UpdateWorldString()
    FTS_Method->UpdateWorldString();
	for(size_t i = 0; i < FTS_Method->_centers.size(); i++){
		EXPECT_NEAR(FTS_Method->_worldstring[mpiid][i], FTS_Method->_centers[i], eps) << "ERROR: _worldstring not updated correctly by UpdateWorldString()";
	}

	// Test StringMethod->GatherNeighbors()
    FTS_Method->SetSendRecvNeighbors(); 
    std::vector<double> lcv0, ucv0;
	lcv0.resize(FTS_Method->_centers.size(), 0);
	ucv0.resize(FTS_Method->_centers.size(), 0);
	// A call to UpdateWorldString() will make the values in worldstring identical to centers
	// Already verified in previous test so this is safe
    FTS_Method->UpdateWorldString(); 
    FTS_Method->GatherNeighbors(&lcv0, &ucv0);
    // Now test that each mpiid has the correct value of centers
    if(FTS_Method->_mpiid == 0)
    {
        for(size_t i = 0; i < FTS_Method->_centers.size(); i++)
        {
            EXPECT_NEAR(lcv0[i], FTS_Method->_worldstring[2][i], eps) << "ERROR: GatherNeighbors() gives incorrect lower CV";
            EXPECT_NEAR(ucv0[i], FTS_Method->_worldstring[1][i], eps) << "ERROR: GatherNeighbors() gives incorrect upper CV";
        }
    }
    else if(FTS_Method->_mpiid == 1)
    {
        for(size_t i = 0; i < FTS_Method->_centers.size(); i++)
        {
            EXPECT_NEAR(lcv0[i], FTS_Method->_worldstring[0][i], eps) << "ERROR: GatherNeighbors() gives incorrect lower CV";
            EXPECT_NEAR(ucv0[i], FTS_Method->_worldstring[2][i], eps) << "ERROR: GatherNeighbors() gives incorrect upper CV";
        }
    }
    if(FTS_Method->_mpiid == 2)
    {
        for(size_t i = 0; i < FTS_Method->_centers.size(); i++)
        {
            EXPECT_NEAR(lcv0[i], FTS_Method->_worldstring[1][i], eps) << "ERROR: GatherNeighbors() gives incorrect lower CV";
            EXPECT_NEAR(ucv0[i], FTS_Method->_worldstring[0][i], eps) << "ERROR: GatherNeighbors() gives incorrect upper CV";
        }
    }

	// Test StringMethod->StringReparam()
    double alphastar = FTS_Method->sqdist(FTS_Method->_centers, lcv0);
    FTS_Method->StringReparam(alphastar); 
    FTS_Method->UpdateWorldString();
    
    //This test is based on empirical results of actually running stringreparam
    if(mpiid == 0)
    {
        EXPECT_NEAR(FTS_Method->_centers[0], -1.1, eps);
        EXPECT_NEAR(FTS_Method->_centers[1], -1.1, eps);
    }
    else if(mpiid == 1)
    {
        EXPECT_NEAR(FTS_Method->_centers[0], 2.38329487, eps);
        EXPECT_NEAR(FTS_Method->_centers[1], -1.6605195328, eps);
    }
    else if(mpiid == 2)
    {
        EXPECT_NEAR(FTS_Method->_centers[0], 2.80, eps);
        EXPECT_NEAR(FTS_Method->_centers[1], 2.2, eps);
    }
	// Test StringMethod->CheckEnd()
	// Check: iteration not met + tol not met
	FTS_Method->_tol = {0.0001, 0.0001};
    FTS_Method->_worldstring = worldstring;
    EXPECT_TRUE(FTS_Method->CheckEnd()) << "ERROR: CheckEnd() exits method earlier than expected";
	// Check: iteration met + tol not met
	FTS_Method->_iteration++;
    EXPECT_DEATH(FTS_Method->CheckEnd(), "") << "ERROR: CheckEnd() fails to exit when max iterations has been met";
	// Check: iteration not met + tol met
	FTS_Method->_iteration--;
	FTS_Method->_tol = {1, 1};
    EXPECT_DEATH(FTS_Method->CheckEnd(), "") << "ERROR: CheckEnd() fails to exit when tolerance criteria has been met";

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


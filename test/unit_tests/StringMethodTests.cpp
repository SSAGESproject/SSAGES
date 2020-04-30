#include "Tests.h"

#include "Methods/FiniteTempString.h"
#include "Methods/Swarm.h"
#include "CVs/MockCV.h"

using namespace SSAGES;

const double sq_dist_eps = 0.001; //Is this too big?  We might not actually get accurate enough results out of square dist
class StringMethodTests : public ::testing::Test
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
				centers_ = {-1.1, -1.1};
			else if(mpiid == 1)
				centers_ = {2.2, -1.8};
			else if(mpiid == 2)
				centers_ = {2.8, 2.2};

			// Set up mock CVs
			CV1 = new MockCV(-0.5, {1, 0, 0}, -2, 2); 
			CV2 = new MockCV(0.5, {0, 1, 0}, -2, 2);
			cvlist.push_back(CV1);
			cvlist.push_back(CV2);
        	
			/* Need a dummy method to test StringMethods functions because we can't create an instance 
			 * of a class with purely virtual functions (as StringMethod has).  
			 * Prefer FTS for consistency when calling common functions. */

        	// Dummy FTS constructor - world, comm, centers, maxiterations, block iterations, tau, cv spring, kappa, spring iter, frequency 
        	FTS_Method = new FiniteTempString(world, comm, centers_, maxiterator, 100, 1, cvspring_, 1, 1, 1);
        	FTS_Method->worldstring_ = worldstring;
        	FTS_Method->mpiid_ = mpiid;
        	FTS_Method->numnodes_ = numnodes;
        	FTS_Method->iteration_ = iteration;
		}

		virtual void TearDown(){
			delete CV1;
			delete CV2;
            
            delete snapshot1;

            delete FTS_Method;
		}

		mxx::comm world;
		mxx::comm comm = world.split(world.rank() < 3 ? world.rank() : 4);
        
        std::vector<double> vec1, vec2, centers_, cvspring_;
		std::vector<std::vector<double>> worldstring;
		unsigned int mpiid, maxiterator, iteration, numnodes;

    	MockCV* CV1;
	    MockCV* CV2;

		CVList cvlist;
        
        Snapshot* snapshot1;
	
        FiniteTempString* FTS_Method;
};

TEST_F(StringMethodTests,dummytest)
{
    //The intention here is to suppress output from functions
    std::cout.setstate(std::ios_base::failbit); 

	// only run this test with 3 ~faux walkers~
	ASSERT_TRUE(3 == world.size()) << "ERROR: String Method unit test must be run with 3 processes" ;

	// Test StringMethod->distance()
	ASSERT_EQ(vec1.size(), vec2.size());
	EXPECT_NEAR(FTS_Method->distance(vec1, vec2), 4.359, sq_dist_eps);
	
	// Test StringMethod->TolCheck()
	// Default is no tolerance defined: tol_ is empty
	EXPECT_FALSE(FTS_Method->TolCheck()) << "ERROR: TolCheck() fails if no tol_ value provided";
	// Initialize tolerance values: begin with tolerance criteria not met 
    FTS_Method->tol_ = {0.0001, 0.0001};
	EXPECT_FALSE(FTS_Method->TolCheck()) << "ERROR: TolCheck() gives false positive";
	// Only one CV within tolerance:
	FTS_Method->tol_[0] = 10;
	EXPECT_FALSE(FTS_Method->TolCheck()) << "ERROR: TolCheck() gives false positive";
	FTS_Method->tol_[1] = 10;
	// After tolerance criteria is met:
	EXPECT_TRUE(FTS_Method->TolCheck()) << "ERROR: TolCheck() gives false negative";
	
	// Test StringMethod->UpdateWorldString()
    FTS_Method->UpdateWorldString(cvlist);
	for(size_t i = 0; i < FTS_Method->centers_.size(); i++){
		EXPECT_NEAR(FTS_Method->worldstring_[mpiid][i], FTS_Method->centers_[i], eps) << "ERROR: worldstring_ not updated correctly by UpdateWorldString()";
	}

	// Test StringMethod->GatherNeighbors()
    FTS_Method->SetSendRecvNeighbors(); 
    std::vector<double> lcv0, ucv0;
	lcv0.resize(FTS_Method->centers_.size(), 0);
	ucv0.resize(FTS_Method->centers_.size(), 0);
	// A call to UpdateWorldString() will make the values in worldstring identical to centers
	// Already verified in previous test so this is safe
    FTS_Method->UpdateWorldString(cvlist); 
    FTS_Method->GatherNeighbors(&lcv0, &ucv0);
    // Now test that each mpiid has the correct value of centers
    if(FTS_Method->mpiid_ == 0)
    {
        for(size_t i = 0; i < FTS_Method->centers_.size(); i++)
        {
            EXPECT_NEAR(lcv0[i], FTS_Method->worldstring_[2][i], eps) << "ERROR: GatherNeighbors() gives incorrect lower CV";
            EXPECT_NEAR(ucv0[i], FTS_Method->worldstring_[1][i], eps) << "ERROR: GatherNeighbors() gives incorrect upper CV";
        }
    }
    else if(FTS_Method->mpiid_ == 1)
    {
        for(size_t i = 0; i < FTS_Method->centers_.size(); i++)
        {
            EXPECT_NEAR(lcv0[i], FTS_Method->worldstring_[0][i], eps) << "ERROR: GatherNeighbors() gives incorrect lower CV";
            EXPECT_NEAR(ucv0[i], FTS_Method->worldstring_[2][i], eps) << "ERROR: GatherNeighbors() gives incorrect upper CV";
        }
    }
    if(FTS_Method->mpiid_ == 2)
    {
        for(size_t i = 0; i < FTS_Method->centers_.size(); i++)
        {
            EXPECT_NEAR(lcv0[i], FTS_Method->worldstring_[1][i], eps) << "ERROR: GatherNeighbors() gives incorrect lower CV";
            EXPECT_NEAR(ucv0[i], FTS_Method->worldstring_[0][i], eps) << "ERROR: GatherNeighbors() gives incorrect upper CV";
        }
    }

	// Test StringMethod->StringReparam()
    double alphastar = FTS_Method->distance(FTS_Method->centers_, lcv0);
    FTS_Method->StringReparam(alphastar); 
    FTS_Method->UpdateWorldString(cvlist);
    
    //This test is based on empirical results of actually running stringreparam
    if(mpiid == 0)
    {
        EXPECT_NEAR(FTS_Method->centers_[0], -1.1, eps);
        EXPECT_NEAR(FTS_Method->centers_[1], -1.1, eps);
    }
    else if(mpiid == 1)
    {
        EXPECT_NEAR(FTS_Method->centers_[0], 2.38329487, eps);
        EXPECT_NEAR(FTS_Method->centers_[1], -1.6605195328, eps);
    }
    else if(mpiid == 2)
    {
        EXPECT_NEAR(FTS_Method->centers_[0], 2.80, eps);
        EXPECT_NEAR(FTS_Method->centers_[1], 2.2, eps);
    }
	// Test StringMethod->CheckEnd()
	// Check: iteration not met + tol not met
	FTS_Method->tol_ = {0.0001, 0.0001};
    FTS_Method->worldstring_ = worldstring;
    EXPECT_TRUE(FTS_Method->CheckEnd(cvlist)) << "ERROR: CheckEnd() exits method earlier than expected";
	// Check: iteration met + tol not met
	FTS_Method->iteration_++;
    EXPECT_DEATH(FTS_Method->CheckEnd(cvlist), "") << "ERROR: CheckEnd() fails to exit when max iterations has been met";
	// Check: iteration not met + tol met
	FTS_Method->iteration_--;
	FTS_Method->tol_ = {1, 1};
    EXPECT_DEATH(FTS_Method->CheckEnd(cvlist), "") << "ERROR: CheckEnd() fails to exit when tolerance criteria has been met";

    std::cout.clear();
}

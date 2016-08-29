#include "../src/Methods/ABF.h"
#include "../src/Snapshot.h"
#include "gtest/gtest.h"
#include <boost/mpi.hpp>

using namespace SSAGES;

// Test calculation up to accuracy of 10^-10
const double eps = 0.0000000001;

class ABFTest : public ::testing::Test {

protected:
    	virtual void SetUp() {
		method = new ABF();
}

	virtual void TearDown() {
		delete method;
}

	ABF *method;

	Snapshot *snapshot;
};

TEST_F(ABFTest,dummytest)
{
}

int main(int argc, char *argv[])
{
    ::testing::InitGoogleTest(&argc, argv);
    int ret = RUN_ALL_TESTS();
    return ret;
}






#include "CVs/MockCV.h"
#include "Snapshot.h"
#include "gtest/gtest.h"
#include <mxx/env.hpp>
#include <mxx/comm.hpp>

using namespace SSAGES;


class MockCVTest : public ::testing::Test {
protected:
    virtual void SetUp() 
    {

        Vector3 grad;
        grad.setZero();

        mocktest = new MockCV(1.2, grad, 5.2, 10.1);
        grad = {1,2,3};
    	mocktest2 = new MockCV(2.3, grad, 6.3, 11.2);
        snapshot1 = new Snapshot(comm, 0);

        auto& pos = snapshot1->GetPositions();
        pos.resize(3);
        auto& ids = snapshot1->GetAtomIDs();
        ids.resize(3);

        for(unsigned int i =0; i <ids.size();i++)
        	ids[i] = i+1;

        pos[0][0] = 0; 
        pos[0][1] = 0; 
        pos[0][2] = 0;
        
        pos[1][0] = 0;
        pos[1][1] = 1;
        pos[1][2] = 0;
        
        pos[2][0] = 1;
        pos[2][1] = 1;
        pos[2][2] = 0;
    }

    virtual void TearDown() {
    	delete mocktest;
    }

    MockCV* mocktest;
    MockCV* mocktest2;

    // Initialize atoms and CV.
    mxx::comm comm;

    Snapshot* snapshot1;

};


TEST_F(MockCVTest, DefaultBehavior)
{
	mocktest->Initialize(*snapshot1);
    mocktest2->Initialize(*snapshot1);

    mocktest->Evaluate(*snapshot1);
	mocktest2->Evaluate(*snapshot1);

	auto& grad1 = mocktest->GetGradient();
	EXPECT_TRUE(grad1[0].isApprox(Vector3{0, 0, 0}));
	EXPECT_TRUE(grad1[1].isApprox(Vector3{0, 0, 0}));
	EXPECT_TRUE(grad1[2].isApprox(Vector3{0, 0, 0}));

	auto& grad2 = mocktest2->GetGradient();
	EXPECT_TRUE(grad2[0].isApprox(Vector3{1, 2, 3}));
	EXPECT_TRUE(grad2[1].isApprox(Vector3{1, 2, 3}));
	EXPECT_TRUE(grad2[2].isApprox(Vector3{1, 2, 3}));

	EXPECT_NEAR(mocktest->GetValue(), 1.2, 1E-8);
	EXPECT_NEAR(mocktest2->GetValue(), 2.3, 1E-8);

    auto& bounds = mocktest->GetBoundaries();
    auto& bounds2 = mocktest2->GetBoundaries();

    EXPECT_NEAR(bounds[0], 5.2, 1E-8);
    EXPECT_NEAR(bounds[1], 10.1, 1E-8);
    EXPECT_NEAR(bounds2[0], 6.3, 1E-8);
    EXPECT_NEAR(bounds2[1], 11.2, 1E-8);
}

int main(int argc, char *argv[])
{
    ::testing::InitGoogleTest(&argc, argv);
    mxx::env env(argc,argv);
    int ret = RUN_ALL_TESTS();
    return ret;
}
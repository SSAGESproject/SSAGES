#include "../src/CVs/TorsionalCV.h"
#include "../src/Snapshot.h"
#include "gtest/gtest.h"
#include <mxx/env.hpp>
#include <mxx/comm.hpp>

using namespace SSAGES;


// 90 degree angle
const double piOver2 = 0.5*M_PI;

class TorsionalCVTest : public ::testing::Test {
protected:
    virtual void SetUp() 
    {

    	tortest = new TorsionalCV(1, 2, 3, 4, true);
    	tortest2 = new TorsionalCV(1, 2, 3, 5, true);
    	tortest3 = new TorsionalCV(1, 2, 3, 6, true);

        snapshot1 = new Snapshot(comm, 0);

        Matrix3 H; 
        H << 100.0, 0.0, 0.0,
             0.0, 100.0, 0.0,
             0.0, 0.0, 100.0;
        snapshot1->SetHMatrix(H);

        snapshot1->SetNumAtoms(6);

        auto& pos = snapshot1->GetPositions();
        pos.resize(6);
        auto& ids = snapshot1->GetAtomIDs();
        ids.resize(6);

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
        
        pos[3][0] = 1;
        pos[3][1] = 2;
        pos[3][2] = 0;

        pos[4][0] = 1;
        pos[4][1] = 1;
        pos[4][2] = 1;

        pos[5][0] = 1;
        pos[5][1] = 0;
        pos[5][2] = 0;

    }

    virtual void TearDown() {
    	delete tortest;
    	delete tortest2;
    	delete tortest3;

    	delete snapshot1;

    }

    TorsionalCV* tortest;
    TorsionalCV* tortest2;
    TorsionalCV* tortest3;

    // Initialize atoms and CV.
    mxx::comm comm;

    Snapshot* snapshot1;

};


TEST_F(TorsionalCVTest, DefaultBehavior)
{
	tortest->Initialize(*snapshot1);
	tortest2->Initialize(*snapshot1);
	tortest3->Initialize(*snapshot1);

	tortest->Evaluate(*snapshot1);
	tortest2->Evaluate(*snapshot1);
	tortest3->Evaluate(*snapshot1);


	auto& grad1 = tortest->GetGradient();
	EXPECT_TRUE(grad1[0].isApprox(Vector3{0, 0, 1}));
	EXPECT_TRUE(grad1[1].isApprox(Vector3{0, 0, -1}));
	EXPECT_TRUE(grad1[2].isApprox(Vector3{0, 0, -1}));
	EXPECT_TRUE(grad1[3].isApprox(Vector3{0, 0, 1}));
	EXPECT_TRUE(grad1[4].isApprox(Vector3{0, 0, 0}));
	EXPECT_TRUE(grad1[5].isApprox(Vector3{0, 0, 0}));

	auto& grad2 = tortest2->GetGradient();
	EXPECT_TRUE(grad2[0].isApprox(Vector3{0, 0, 1}));
	EXPECT_TRUE(grad2[1].isApprox(Vector3{0, 0, -1}));
	EXPECT_TRUE(grad2[2].isApprox(Vector3{0, 1, 0}));
	EXPECT_TRUE(grad2[3].isApprox(Vector3{0, 0, 0}));
	EXPECT_TRUE(grad2[4].isApprox(Vector3{0, -1, 0}));
	EXPECT_TRUE(grad2[5].isApprox(Vector3{0, 0, 0}));
	
	auto& grad3 = tortest3->GetGradient();
	EXPECT_TRUE(grad3[0].isApprox(Vector3{0, 0, 1}));
	EXPECT_TRUE(grad3[1].isApprox(Vector3{0, 0, -1}));
	EXPECT_TRUE(grad3[2].isApprox(Vector3{0, 0, 1}));
	EXPECT_TRUE(grad3[3].isApprox(Vector3{0, 0, 0}));
	EXPECT_TRUE(grad3[4].isApprox(Vector3{0, 0, 0}));
	EXPECT_TRUE(grad3[5].isApprox(Vector3{0, 0, -1}));

	EXPECT_NEAR(tortest->GetValue(), 3.14159, 0.01);
	EXPECT_NEAR(tortest2->GetValue(), -1.570796, 0.01);
	EXPECT_NEAR(tortest3->GetValue(), 0, 0.01);
}

int main(int argc, char *argv[])
{
    ::testing::InitGoogleTest(&argc, argv);
    mxx::env env(argc,argv);
    int ret = RUN_ALL_TESTS();
    return ret;
}

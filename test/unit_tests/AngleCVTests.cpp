#include "../src/CVs/AngleCV.h"
#include "../src/Snapshot.h"
#include "gtest/gtest.h"
#include <mxx/env.hpp>
#include <mxx/comm.hpp>

using namespace SSAGES;

// Test up to accuracy of 10^-10
const double eps = 1e-10;

// 90 degree angle
const double piOver2 = 0.5*M_PI;

class AngleCVTest : public ::testing::Test {
protected:
    virtual void SetUp() 
    {

        angletest = new AngleCV(1, 2, 3);
        angletest2 = new AngleCV(1, 3, 4);
        angletest3 = new AngleCV(3, 2, 4);
        angletest4 = new AngleCV(2, 1, 5);

        snapshot1 = new Snapshot(comm, 0);

        auto& pos = snapshot1->GetPositions();
        pos.resize(5);
        auto& ids = snapshot1->GetAtomIDs();
        ids.resize(5);

        Matrix3 H; 
        H << 100.0, 0.0, 0.0,
             0.0, 100.0, 0.0,
             0.0, 0.0, 100.0;

        snapshot1->SetHMatrix(H);
        snapshot1->SetNumAtoms(5);

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

        pos[4][0] = 0;
        pos[4][1] = -1;
        pos[4][2] = 0;

    }

    virtual void TearDown() {
    	delete angletest;
    	delete angletest2;
    	delete angletest3;
    	delete angletest4;

    	delete snapshot1;

    }

    AngleCV* angletest;
    AngleCV* angletest2;
    AngleCV* angletest3;
    AngleCV* angletest4;

    // Initialize atoms and CV.
    mxx::comm comm;

    Snapshot* snapshot1;

};

TEST_F(AngleCVTest, DefaultBehavior)
{
    angletest->Initialize(*snapshot1);
    angletest2->Initialize(*snapshot1);
    angletest3->Initialize(*snapshot1);
    angletest4->Initialize(*snapshot1);

    angletest->Evaluate(*snapshot1);
    angletest2->Evaluate(*snapshot1);
    angletest3->Evaluate(*snapshot1);
    angletest4->Evaluate(*snapshot1);

    EXPECT_NEAR(angletest->GetValue(), 1.570796, 0.01);
    EXPECT_NEAR(angletest2->GetValue(), 2.356, 0.01);
    EXPECT_NEAR(angletest3->GetValue(), 0.785398, 0.01);
    EXPECT_NEAR(angletest4->GetValue(), 3.14159, 0.01);
}


int main(int argc, char *argv[])
{
    ::testing::InitGoogleTest(&argc, argv);
    mxx::env env(argc,argv);
    int ret = RUN_ALL_TESTS();
    return ret;
}
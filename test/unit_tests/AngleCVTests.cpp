#include "Tests.h"

#include "CVs/AngleCV.h"

using namespace SSAGES;

class AngleCVTests : public ::testing::Test {
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

TEST_F(AngleCVTests, DefaultBehavior)
{
    angletest->Initialize(*snapshot1);
    angletest2->Initialize(*snapshot1);
    angletest3->Initialize(*snapshot1);
    angletest4->Initialize(*snapshot1);

    angletest->Evaluate(*snapshot1);
    angletest2->Evaluate(*snapshot1);
    angletest3->Evaluate(*snapshot1);
    angletest4->Evaluate(*snapshot1);

    EXPECT_NEAR(angletest->GetValue(), M_PI/2, eps);
    EXPECT_NEAR(angletest2->GetValue(), 3*M_PI/4, eps);
    EXPECT_NEAR(angletest3->GetValue(), M_PI/4, eps);
    EXPECT_NEAR(angletest4->GetValue(), M_PI, eps);
}

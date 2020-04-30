#include "Tests.h"

#include "CVs/RouseModeCV.h"

using namespace SSAGES;

class RouseModeCVTests : public ::testing::Test { 
protected:
	virtual void SetUp()
	{
        X0 = new RouseModeCV({{1},{2},{3},{4}},0); // 0th Rouse mode
        X1 = new RouseModeCV({{1},{2},{3},{4}},1); // 1st Rouse mode

    	// Set up snapshot 
        snapshot1 = new Snapshot(comm, 0);
        snapshot1->SetNumAtoms(4);
        Matrix3 H;
        H << 10, 0, 0, 0, 10, 0, 0, 0, 10;
        snapshot1->SetHMatrix(H);
        

        // Atom positions
        auto& pos1 = snapshot1->GetPositions();
        pos1.resize(4);
        pos1[0][0] = 1.0;
        pos1[0][1] = 0.0;
        pos1[0][2] = 0.0;

        pos1[1][0] = 2.0;
        pos1[1][1] = 0.0;
        pos1[1][2] = 0.0;

        pos1[2][0] = 3.0;
        pos1[2][1] = 0.0;
        pos1[2][2] = 0.0;

        pos1[3][0] = 4.0;
        pos1[3][1] = 0.0;
        pos1[3][2] = 0.0;

        // Atom IDs
        auto& ids1 = snapshot1->GetAtomIDs();
        ids1.resize(4);
        ids1[0] = 1;
        ids1[1] = 2;
        ids1[2] = 3;
        ids1[3] = 4;
 
        // Atom Masses
        auto& mass = snapshot1->GetMasses();
        mass.resize(4);
        mass[0] = 1;
        mass[1] = 1;
        mass[2] = 1;
        mass[3] = 1;
	}

	virtual void TearDown() 
	{
		delete X0; 
		delete X1; 
		delete snapshot1;
	}


    RouseModeCV* X0;
    RouseModeCV* X1;
	mxx::comm comm;
	Snapshot* snapshot1;
};

TEST_F(RouseModeCVTests, DefaultBehavior)
{
    X0->Initialize(*snapshot1);
    X1->Initialize(*snapshot1);

    X0->Evaluate(*snapshot1);
    X1->Evaluate(*snapshot1);

    // CHECK VALUES
    EXPECT_NEAR(X0->GetValue(),5.00000    , eps); 
    EXPECT_NEAR(X1->GetValue(),2.23044250 , eps);

    // CHECK GRADIENTS
    // 0th mode
    EXPECT_NEAR(X0->GetGradient()[0][0],0.500, eps);
    EXPECT_NEAR(X0->GetGradient()[0][1],0.000 , eps);
    EXPECT_NEAR(X0->GetGradient()[0][2],0.000 , eps);

    EXPECT_NEAR(X0->GetGradient()[1][0],0.500 , eps);
    EXPECT_NEAR(X0->GetGradient()[1][1],0.000 , eps);
    EXPECT_NEAR(X0->GetGradient()[1][2],0.000 , eps);
    
    EXPECT_NEAR(X0->GetGradient()[2][0],0.500 , eps);
    EXPECT_NEAR(X0->GetGradient()[2][1],0.000 , eps);
    EXPECT_NEAR(X0->GetGradient()[2][2],0.000 , eps);

    EXPECT_NEAR(X0->GetGradient()[3][0],0.500 , eps);
    EXPECT_NEAR(X0->GetGradient()[3][1],0.000 , eps);
    EXPECT_NEAR(X0->GetGradient()[3][2],0.000 , eps);

    // 1st mode
    EXPECT_NEAR(X1->GetGradient()[0][0],-0.65328148, eps);
    EXPECT_NEAR(X1->GetGradient()[0][1],0.000 , eps);
    EXPECT_NEAR(X1->GetGradient()[0][2],0.000 , eps);

    EXPECT_NEAR(X1->GetGradient()[1][0],-0.27059805 , eps);
    EXPECT_NEAR(X1->GetGradient()[1][1],0.000 , eps);
    EXPECT_NEAR(X1->GetGradient()[1][2],0.000 , eps);
    
    EXPECT_NEAR(X1->GetGradient()[2][0],0.27059805 , eps);
    EXPECT_NEAR(X1->GetGradient()[2][1],0.000 , eps);
    EXPECT_NEAR(X1->GetGradient()[2][2],0.000 , eps);

    EXPECT_NEAR(X1->GetGradient()[3][0],0.65328148 , eps);
    EXPECT_NEAR(X1->GetGradient()[3][1],0.000 , eps);
    EXPECT_NEAR(X1->GetGradient()[3][2],0.000 , eps);

}

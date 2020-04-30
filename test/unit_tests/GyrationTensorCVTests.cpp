#include "Tests.h"

#include "CVs/GyrationTensorCV.h"

using namespace SSAGES; 

class GyrationTensorCVTests : public ::testing::Test { 
protected:
    virtual void SetUp()
    {
        cv = new GyrationTensorCV({1,2,3}, Rg);
        cv_xy = new GyrationTensorCV({1,2,3}, Rg, true, true, false);

        // Set up snapshot No. 1
        snapshot1 = new Snapshot(comm, 0);
        
        Matrix3 H;
        H << 100, 0, 0, 0, 100, 0, 0, 0, 100;
        snapshot1->SetHMatrix(H);
        
        snapshot1->SetNumAtoms(3);

        auto& pos1 = snapshot1->GetPositions();
        pos1.resize(3);
        pos1[0][0] = 0.5;
        pos1[0][1] = 0.5;
        pos1[0][2] = 0.5;

        pos1[1][0] = 1.0;
        pos1[1][1] = 1.0;
        pos1[1][2] = 1.0;

        pos1[2][0] = 1.5;
        pos1[2][1] = 1.5;
        pos1[2][2] = 1.5;

        auto& ids1 = snapshot1->GetAtomIDs();
        ids1.resize(3);
        ids1[0] = 1;
        ids1[1] = 2;
        ids1[2] = 3;
 
        auto& mass = snapshot1->GetMasses();
        mass.resize(3);
        mass[0] = 1;
        mass[1] = 1;
        mass[2] = 1;
    }

    virtual void TearDown() 
    {
        delete cv;
        delete cv_xy;
        delete snapshot1;
    }

    GyrationTensorCV* cv;
    GyrationTensorCV* cv_xy;
    mxx::comm comm;
    Snapshot* snapshot1;
};

TEST_F(GyrationTensorCVTests, CompareRg)
{
    cv->Initialize(*snapshot1);
    cv->Evaluate(*snapshot1);
    EXPECT_NEAR(cv->GetValue(), 0.5, eps);
    EXPECT_NEAR(cv->GetGradient()[0][0], -0.2222222222, eps);
    EXPECT_NEAR(cv->GetGradient()[0][1], -0.2222222222, eps);
    EXPECT_NEAR(cv->GetGradient()[0][2], -0.2222222222, eps);

    EXPECT_NEAR(cv->GetGradient()[1][0], 0., eps);
    EXPECT_NEAR(cv->GetGradient()[1][1], 0., eps);
    EXPECT_NEAR(cv->GetGradient()[1][2], 0., eps);

    EXPECT_NEAR(cv->GetGradient()[2][0], 0.2222222222, eps);
    EXPECT_NEAR(cv->GetGradient()[2][1], 0.2222222222, eps);
    EXPECT_NEAR(cv->GetGradient()[2][2], 0.2222222222, eps);

    cv_xy->Initialize(*snapshot1);
    cv_xy->Evaluate(*snapshot1);
    EXPECT_NEAR(cv_xy->GetValue(), 0.3333333333, eps);
    EXPECT_NEAR(cv_xy->GetGradient()[0][0], -0.2222222222, eps);
    EXPECT_NEAR(cv_xy->GetGradient()[0][1], -0.2222222222, eps);
    EXPECT_NEAR(cv_xy->GetGradient()[0][2], 0., eps);

    EXPECT_NEAR(cv_xy->GetGradient()[1][0], 0., eps);
    EXPECT_NEAR(cv_xy->GetGradient()[1][1], 0., eps);
    EXPECT_NEAR(cv_xy->GetGradient()[1][2], 0., eps);

    EXPECT_NEAR(cv_xy->GetGradient()[2][0], 0.2222222222, eps);
    EXPECT_NEAR(cv_xy->GetGradient()[2][1], 0.2222222222, eps);
    EXPECT_NEAR(cv_xy->GetGradient()[2][2], 0., eps);
}

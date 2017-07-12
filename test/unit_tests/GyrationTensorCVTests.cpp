#include "CVs/GyrationTensorCV.h"
#include "Snapshot.h"
#include "gtest/gtest.h"
#include <mxx/env.hpp>
#include <mxx/comm.hpp>

using namespace SSAGES; 

// Test accuracy up to 1e-10. 
constexpr double eps = 1e-10;

class GyrationTensorCVTests : public ::testing::Test { 
protected:
	virtual void SetUp()
	{
    	cv = new GyrationTensorCV({1,2,3}, Rg);

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
		delete snapshot1;
	}

	GyrationTensorCV* cv; 
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
}

int main(int argc, char *argv[])
{
    ::testing::InitGoogleTest(&argc, argv);
    mxx::env env(argc,argv);
    int ret = RUN_ALL_TESTS();
    return ret;
}
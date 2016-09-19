#include "../src/CVs/RadiusOfGyrationCV.h"
#include "../src/Snapshot.h"
#include <iostream>
#include <fstream>
#include "gtest/gtest.h"
#include <boost/mpi.hpp>

using namespace SSAGES;

// Test up to accuracy of 10^-10
const double eps = 1e-10;

class RgCVTests : public ::testing::Test {
protected:
    virtual void SetUp() 
    {
    	Rg = new RadiusOfGyrationCV({1,2,3});

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

    RadiusOfGyrationCV *Rg;
    boost::mpi::communicator comm;

    Snapshot *snapshot1;
};

TEST_F(RgCVTests, CompareRg)
{
	Rg->Initialize(*snapshot1);
    Rg->Evaluate(*snapshot1);
    EXPECT_NEAR(Rg->GetValue(), 0.70710678118, eps);
}


int main(int argc, char *argv[])
{
    ::testing::InitGoogleTest(&argc, argv);
    boost::mpi::environment env(argc,argv);
    int ret = RUN_ALL_TESTS();
    return ret;
}
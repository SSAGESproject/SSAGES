#include "../src/CVs/ParticleSeparationCV.h"
#include "../src/Snapshot.h"
#include <iostream>
#include <fstream>
#include "gtest/gtest.h"
#include <boost/mpi.hpp>

using namespace SSAGES;

// Test up to accuracy of 10^-10
const double eps = 1e-10;

class ParticleSeparationCVTests : public ::testing::Test {
protected:
    virtual void SetUp() 
    {

        COMDist = new ParticleSeparationCV({1,2,3}, {4,5,6});

    	// Set up snapshot No. 1
        // Snapshot 1 contains six atoms
        snapshot1 = new Snapshot(comm, 0);

        Matrix3 H; 
        H << 100.0, 0.0, 0.0,
             0.0, 100.0, 0.0,
             0.0, 0.0, 100.0;
        snapshot1->SetHMatrix(H);
        snapshot1->SetNumAtoms(6);

        auto& map = snapshot1->GetIDMap();
        map[1] = 0;
        map[2] = 1;
        map[3] = 2;
        map[4] = 3;
        map[5] = 4;
        map[6] = 5;

        auto& pos1 = snapshot1->GetPositions();
        pos1.resize(6);
        pos1[0][0] = 0.5;
        pos1[0][1] = 1.5;
        pos1[0][2] = 0.5;

        pos1[1][0] = 1.0;
        pos1[1][1] = 0.5;
        pos1[1][2] = 2.0;

        pos1[2][0] = 1.5;
        pos1[2][1] = 1.5;
        pos1[2][2] = 1.5;

        pos1[3][0] = 2.5;
        pos1[3][1] = 4.5;
        pos1[3][2] = 6.5;

        pos1[4][0] = 3.0;
        pos1[4][1] = 8.0;
        pos1[4][2] = 6.0;

        pos1[5][0] = 3.5;
        pos1[5][1] = 4.5;
        pos1[5][2] = 5.5;

        auto& ids1 = snapshot1->GetAtomIDs();
        ids1.resize(6);
        ids1[0] = 1;
        ids1[1] = 2;
        ids1[2] = 3;
        ids1[3] = 4;
        ids1[4] = 5;
        ids1[5] = 6;

        auto& mass = snapshot1->GetMasses();
        mass.resize(6);
        mass[0] = 1;
        mass[1] = 2;
        mass[2] = 3;
        mass[3] = 4;
        mass[4] = 5;
        mass[5] = 6;

    }

    ParticleSeparationCV *COMDist;

    boost::mpi::communicator comm;

    Snapshot *snapshot1;
};

TEST_F(ParticleSeparationCVTests, compareCOMDist)
{
	COMDist->Initialize(*snapshot1);
    COMDist->Evaluate(*snapshot1);

    EXPECT_NEAR(COMDist->GetValue(), 6.596547918756026, eps);
}

int main(int argc, char *argv[])
{
    ::testing::InitGoogleTest(&argc, argv);
    boost::mpi::environment env(argc,argv);
    int ret = RUN_ALL_TESTS();
    return ret;
}
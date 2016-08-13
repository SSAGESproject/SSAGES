#include "../src/CVs/RMSDCV.h"
#include "../src/Snapshot.h"
#include <iostream>
#include <fstream>
#include "gtest/gtest.h"
#include <boost/mpi.hpp>

using namespace SSAGES;

// Test up to accuracy of 10^-10
const double eps = 1e-10;

// 90 degree angle
const double piOver2 = 0.5*M_PI;

class RMSDCVTest : public ::testing::Test {
protected:
    virtual void SetUp() 
    {
        //Set up xyz files
        std::ofstream myfile;

        //Normal file with extra line (should pass)
        myfile.open (filexyz1);
        myfile << "3\n";
        myfile << "Comments here\n";
        myfile << "1 0.0 2.2 0.8\n";
        myfile << "2 1.3 2.5 0.5\n";
        myfile << "3 1.3 1.2 1.3";
        myfile.close();

    	normal_RMSD = new RMSDCV({1,3}, filexyz1, true);
    	no_range_RMSD = new RMSDCV({1,2,3}, filexyz1, false);

    	// Set up snapshot No. 1
        // Snapshot 1 contains three atoms
        snapshot1 = new Snapshot(comm, 0);

        auto& latconst = snapshot1->GetLatticeConstants();

        latconst[0] = 100.0;
        latconst[1] = 100.0;
        latconst[2] = 100.0;
        latconst[3] = piOver2;
        latconst[4] = piOver2;
        latconst[5] = piOver2;

        auto& pos1 = snapshot1->GetPositions();
        pos1.resize(3);
        pos1[0][0] = 0.0;
        pos1[0][1] = 2.2;
        pos1[0][2] = 0.9;

        pos1[1][0] = 1.3;
        pos1[1][1] = 2.1;
        pos1[1][2] = 0.4;

        pos1[2][0] = 1.3;
        pos1[2][1] = 1.0;
        pos1[2][2] = 1.1;

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

    virtual void TearDown() {

        std::remove(filexyz1.c_str());
    	delete normal_RMSD;
    	delete no_range_RMSD;

    	delete snapshot1;
    }

    std::string filexyz1 = "test1.xyz";
    RMSDCV *normal_RMSD;
    RMSDCV *no_range_RMSD;

    boost::mpi::communicator comm;

    Snapshot *snapshot1;
};

TEST_F(RMSDCVTest, compareRMSD)
{
	normal_RMSD->Initialize(*snapshot1);


    EXPECT_NEAR(normal_RMSD->GetValue(), 0.29439202887, 1E-8);
}

TEST_F(RMSDCVTest, TestRange)
{
	normal_RMSD->Initialize(*snapshot1);
	no_range_RMSD->Initialize(*snapshot1);

	EXPECT_EQ(normal_RMSD->GetValue(), no_range_RMSD->GetValue());
}

TEST_F(RMSDCVTest, BadRange)
{
	EXPECT_ANY_THROW(new RMSDCV({1,2,3}, filexyz1, true));
}

int main(int argc, char *argv[])
{
    ::testing::InitGoogleTest(&argc, argv);
    int ret = RUN_ALL_TESTS();
    return ret;
}
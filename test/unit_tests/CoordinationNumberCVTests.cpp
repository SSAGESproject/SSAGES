#include "../src/CVs/CoordinationNumberCV.h"
#include "../src/Snapshot.h"
#include "gtest/gtest.h"
#include <boost/mpi.hpp>

using namespace SSAGES;

// Test up to accuracy of 10^-10
constexpr double eps = 1e-10;

class CoordinationNumberCVTest : public ::testing::Test{
protected:
    virtual void SetUp()
    {
        cv1 = new CoordinationNumberCV({1, 2, 3}, {4}, SwitchingFunction(0, 1.5, 12, 24));
        cv2 = new CoordinationNumberCV({4}, {1, 2, 3}, SwitchingFunction(0, 1.5, 12, 24));

        snapshot = new Snapshot(comm, 0);
        Matrix3 H; 
        H << 100.0, 0.0, 0.0,
             0.0, 100.0, 0.0,
             0.0, 0.0, 100.0;
        snapshot->SetHMatrix(H);
        snapshot->SetNumAtoms(4);

        auto& ids = snapshot->GetAtomIDs();
        auto& pos = snapshot->GetPositions();
        
        if(comm.size() == 2)
        {
            pos.resize(2);
            ids.resize(2);
            if(comm.rank() == 0)
            {
                pos[0][0] = 0.5;
                pos[0][1] = 1.5;
                pos[0][2] = 0.5;

                pos[1][0] = 1.0;
                pos[1][1] = 0.5;
                pos[1][2] = 2.0;
                
                ids[0] = 1;
                ids[1] = 2;
            }
            else
            {
                pos[0][0] = 1.5;
                pos[0][1] = 1.5;
                pos[0][2] = 1.5;

                pos[1][0] = 1.5;
                pos[1][1] = 0.5;
                pos[1][2] = 1.5;   

                ids[0] = 3;
                ids[1] = 4;
            }
        }
        else
        {
            pos.resize(4);
            pos[0][0] = 0.5;
            pos[0][1] = 1.5;
            pos[0][2] = 0.5;

            pos[1][0] = 1.0;
            pos[1][1] = 0.5;
            pos[1][2] = 2.0;

            pos[2][0] = 1.5;
            pos[2][1] = 1.5;
            pos[2][2] = 1.5;

            pos[3][0] = 1.5;
            pos[3][1] = 0.5;
            pos[3][2] = 1.5;
            
            ids.resize(4);
            ids[0] = 1;
            ids[1] = 2;
            ids[2] = 3;
            ids[3] = 4;
        }
    }

    virtual void TearDown()
    {
        delete snapshot; 
        delete cv1;
        delete cv2;
    }

    Snapshot* snapshot; 
    CoordinationNumberCV *cv1, *cv2;
    boost::mpi::communicator comm;
};

TEST_F(CoordinationNumberCVTest, DefaultBehavior)
{
    cv1->Initialize(*snapshot);
    cv1->Evaluate(*snapshot);
    EXPECT_NEAR(cv1->GetValue(), 2.143319272335021, eps);

    cv2->Initialize(*snapshot);
    cv2->Evaluate(*snapshot);
    EXPECT_NEAR(cv2->GetValue(), 2.143319272335021, eps);
}

int main(int argc, char *argv[])
{
    ::testing::InitGoogleTest(&argc, argv);
    boost::mpi::environment env(argc,argv);
    int ret = RUN_ALL_TESTS();
    return ret;
}

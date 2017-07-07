#include "../src/CVs/CoordinationNumberCV.h"
#include "../src/Snapshot.h"
#include "gtest/gtest.h"
#include <mxx/env.hpp>
#include <mxx/comm.hpp>

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
    mxx::comm comm;
};

TEST_F(CoordinationNumberCVTest, DefaultBehavior)
{
    cv1->Initialize(*snapshot);
    cv1->Evaluate(*snapshot);
    EXPECT_NEAR(cv1->GetValue(), 2.143319272335021, eps);

    auto& grad1 = cv1->GetGradient(); 
    // NOTE: Gradient checks are only valid for single processor
    // but value checks for 2. This needs to be adjusted later. 
    if(comm.size() == 1)
    {
        EXPECT_NEAR(grad1[0][0], 0.5130418964267503, eps);
        EXPECT_NEAR(grad1[0][1], -0.5130418964267503, eps);
        EXPECT_NEAR(grad1[0][2], 0.5130418964267503, eps);

        EXPECT_NEAR(grad1[1][0], 0.0014447794902723533, eps);
        EXPECT_NEAR(grad1[1][1], 0, eps);
        EXPECT_NEAR(grad1[1][2], -0.0014447794902723533, eps);

        EXPECT_NEAR(grad1[2][0], 0, eps);
        EXPECT_NEAR(grad1[2][1], -0.09107879745469, eps);
        EXPECT_NEAR(grad1[2][2], 0, eps);

        EXPECT_NEAR(grad1[3][0], -0.51448667591702268, eps);
        EXPECT_NEAR(grad1[3][1], 0.60412069388144074, eps);
        EXPECT_NEAR(grad1[3][2], -0.51159711693647802, eps);
    }

    // Gradient should sum to zero.
    double sum = 0.; 
    for(auto& v : grad1)
        sum += v.array().sum();
    sum = mxx::allreduce(sum, std::plus<double>(), comm);
    EXPECT_NEAR(sum, 0, eps);

    // CV2 should match CV1.
    cv2->Initialize(*snapshot);
    cv2->Evaluate(*snapshot);

    EXPECT_NEAR(cv2->GetValue(), cv1->GetValue(), eps);

    auto& grad2 = cv2->GetGradient();
    sum = 0.; 
    for(auto& v : grad2)
        sum += v.array().sum();
    sum = mxx::allreduce(sum, std::plus<double>(), comm);
    EXPECT_NEAR(sum, 0, eps);

    if(comm.size() == 1)
    {
        EXPECT_NEAR(grad2[0][0], grad1[0][0], eps);
        EXPECT_NEAR(grad2[0][1], grad1[0][1], eps);
        EXPECT_NEAR(grad2[0][2], grad1[0][2], eps);

        EXPECT_NEAR(grad2[1][0], grad1[1][0], eps);
        EXPECT_NEAR(grad2[1][1], grad1[1][1], eps);
        EXPECT_NEAR(grad2[1][2], grad1[1][2], eps);

        EXPECT_NEAR(grad2[2][0], grad1[2][0], eps);
        EXPECT_NEAR(grad2[2][1], grad1[2][1], eps);
        EXPECT_NEAR(grad2[2][2], grad1[2][2], eps);

        EXPECT_NEAR(grad2[3][0], grad1[3][0], eps);
        EXPECT_NEAR(grad2[3][1], grad1[3][1], eps);
        EXPECT_NEAR(grad2[3][2], grad1[3][2], eps);
    }
}

int main(int argc, char *argv[])
{
    ::testing::InitGoogleTest(&argc, argv);
    mxx::env env(argc,argv);
    int ret = RUN_ALL_TESTS();
    return ret;
}

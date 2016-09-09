#include "../src/CVs/ParticlePositionCV.h"
#include "../src/Snapshot.h"
#include "gtest/gtest.h"
#include <boost/mpi.hpp>

using namespace SSAGES;

// Test calculation up to accuracy of 10^-10
const double eps = 0.0000000001;

class ParticlePositionCVTest : public ::testing::Test {

protected:
    virtual void SetUp() {
        // Set up various Atom position CVs
        fullyFunctional = new ParticlePositionCV({1}, Vector3({1.0, 1.0, 1.0}), true, true, true);
        constrainedInX = new ParticlePositionCV({1}, Vector3({3.0, 2.1, 4.6}), false, true, true);
        fullyConstrained = new ParticlePositionCV({1}, Vector3({1.0, 1.0, 1.0}), false, false, false);
        illegalAtomID = new ParticlePositionCV({-3}, Vector3({1.0, 1.0, 1.0}), true, true, true);

        // Set up snapshot No. 1
        // Snapshot 1 contains two atoms
        snapshot1 = new Snapshot(comm, 0);
        Matrix3 H;
        H << 100, 0, 0, 0, 100, 0, 0, 0, 100;
        snapshot1->SetHMatrix(H);

        auto& pos1 = snapshot1->GetPositions();
        pos1.resize(2);
        pos1[0][0] = -1.2;
        pos1[0][1] = 2.2;
        pos1[0][2] = 0.9;

        pos1[1][0] = 1.3;
        pos1[1][1] = 2.1;
        pos1[1][2] = 0.4;

        auto& ids1 = snapshot1->GetAtomIDs();
        ids1.resize(2);
        ids1[0] = 1;
        ids1[1] = 2;

        auto& map1 = snapshot1->GetIDMap();
        map1[1] = 0;
        map1[2] = 1; 

        auto& mass1 = snapshot1->GetMasses();
        mass1.push_back(1.0);

        snapshot1->SetNumAtoms(2);

        // Set up snapshot No. 2
        // Positions of snapshot 2 are identical to CV position
        snapshot2 = new Snapshot(comm, 0);
        snapshot2->SetHMatrix(H);

        auto& map2 = snapshot2->GetIDMap();
        map2[1] = 0;

        auto& pos2 = snapshot2->GetPositions();
        pos2.resize(1);
        pos2[0][0] = 1.0;
        pos2[0][1] = 1.0;
        pos2[0][2] = 1.0;

        auto& ids2 = snapshot2->GetAtomIDs();
        ids2.resize(1);
        ids2[0] = 1;

        auto& mass2 = snapshot2->GetMasses();
        mass2.push_back(1.0);
        snapshot2->SetNumAtoms(1);
    }

    virtual void TearDown() {
        delete fullyFunctional;
        delete constrainedInX;
        delete fullyConstrained;
        delete illegalAtomID;

        delete snapshot1;
        delete snapshot2;
    }

    ParticlePositionCV *fullyFunctional;
    ParticlePositionCV *constrainedInX;
    ParticlePositionCV *fullyConstrained;
    ParticlePositionCV *illegalAtomID;

    boost::mpi::communicator comm;

    Snapshot *snapshot1;
    Snapshot *snapshot2;
};

TEST_F(ParticlePositionCVTest, ValueTest1) {
    // Initialize CVs
    ASSERT_NO_THROW(fullyFunctional->Initialize(*snapshot1));
    ASSERT_NO_THROW(constrainedInX->Initialize(*snapshot1));
    ASSERT_NO_THROW(fullyConstrained->Initialize(*snapshot1));
    ASSERT_ANY_THROW(illegalAtomID->Initialize(*snapshot1));

    // Evaluate CVs
    fullyFunctional->Evaluate(*snapshot1);
    constrainedInX->Evaluate(*snapshot1);
    fullyConstrained->Evaluate(*snapshot1);

    // Test if CV values are as expected
    EXPECT_NEAR(fullyFunctional->GetValue(), 2.5079872408, eps);
    EXPECT_NEAR(constrainedInX->GetValue(), 3.7013511047, eps);
    EXPECT_DOUBLE_EQ(fullyConstrained->GetValue(), 0.0);
}

TEST_F(ParticlePositionCVTest, ValueTest2) {
    // Initialize CVs
    fullyFunctional->Initialize(*snapshot2);
    constrainedInX->Initialize(*snapshot2);
    fullyConstrained->Initialize(*snapshot2);

    // Evaluate CVs
    fullyFunctional->Evaluate(*snapshot2);
    constrainedInX->Evaluate(*snapshot2);
    fullyConstrained->Evaluate(*snapshot2);

    // Test if CV values are as expected
    EXPECT_DOUBLE_EQ(fullyFunctional->GetValue(), 0.0);
    EXPECT_NEAR(constrainedInX->GetValue(), 3.7643060449, eps);
    EXPECT_DOUBLE_EQ(fullyConstrained->GetValue(), 0.0);
}

TEST_F(ParticlePositionCVTest, PeriodicValueTest) {
    // Set test values
    double value1 = -15.3;
    double value2 = 0.0;
    double value3 = 4.2;
    double value4 = 1928437.3;

    EXPECT_DOUBLE_EQ(fullyFunctional->GetPeriodicValue(value1), value1);
    EXPECT_DOUBLE_EQ(fullyFunctional->GetPeriodicValue(value2), value2);
    EXPECT_DOUBLE_EQ(fullyFunctional->GetPeriodicValue(value3), value3);
    EXPECT_DOUBLE_EQ(fullyFunctional->GetPeriodicValue(value4), value4);
}

TEST_F(ParticlePositionCVTest, GradientTest1) {
    // Initialize CVs
    fullyFunctional->Initialize(*snapshot1);
    constrainedInX->Initialize(*snapshot1);
    fullyConstrained->Initialize(*snapshot1);

    // Evaluate CVs
    fullyFunctional->Evaluate(*snapshot1);
    constrainedInX->Evaluate(*snapshot1);
    fullyConstrained->Evaluate(*snapshot1);

    // Test if CV values are as expected

    // Test fully functional CV
    std::vector<Vector3> fullyFunctionalGradient = fullyFunctional->GetGradient();
    Vector3 gradient1 = fullyFunctionalGradient.at(0);
    Vector3 gradient2 = fullyFunctionalGradient.at(1);

    EXPECT_NEAR(gradient1[0], -0.8771974451, eps);
    EXPECT_NEAR(gradient1[1], 0.4784713337, eps);
    EXPECT_NEAR(gradient1[2], -0.0398726111, eps);

    EXPECT_DOUBLE_EQ(gradient2[0], 0.0);
    EXPECT_DOUBLE_EQ(gradient2[1], 0.0);
    EXPECT_DOUBLE_EQ(gradient2[2], 0.0);

    // Test CV constrained in x dimension
    std::vector<Vector3> constrainedInXGradient = constrainedInX->GetGradient();
    gradient1 = constrainedInXGradient.at(0);
    gradient2 = constrainedInXGradient.at(1);

    EXPECT_DOUBLE_EQ(gradient1[0], 0.0);
    EXPECT_NEAR(gradient1[1], 0.0270171613, eps);
    EXPECT_NEAR(gradient1[2], -0.9996349699, eps);

    EXPECT_DOUBLE_EQ(gradient2[0], 0.0);
    EXPECT_DOUBLE_EQ(gradient2[1], 0.0);
    EXPECT_DOUBLE_EQ(gradient2[2], 0.0);

    // Test fully constrained CV
    std::vector<Vector3> fullyConstrainedGradient = fullyConstrained->GetGradient();
    gradient1 = fullyConstrainedGradient.at(0);
    gradient2 = fullyConstrainedGradient.at(1);

    EXPECT_DOUBLE_EQ(gradient1[0], 0.0);
    EXPECT_DOUBLE_EQ(gradient1[1], 0.0);
    EXPECT_DOUBLE_EQ(gradient1[2], 0.0);

    EXPECT_DOUBLE_EQ(gradient2[0], 0.0);
    EXPECT_DOUBLE_EQ(gradient2[1], 0.0);
    EXPECT_DOUBLE_EQ(gradient2[2], 0.0);
}

TEST_F(ParticlePositionCVTest, GradientTest2) {
    // Initialize CVs
    fullyFunctional->Initialize(*snapshot2);
    constrainedInX->Initialize(*snapshot2);
    fullyConstrained->Initialize(*snapshot2);

    // Evaluate CVs
    fullyFunctional->Evaluate(*snapshot2);
    constrainedInX->Evaluate(*snapshot2);
    fullyConstrained->Evaluate(*snapshot2);

    // Test if CV values are as expected

    // Test fully functional CV
    std::vector<Vector3> fullyFunctionalGradient = fullyFunctional->GetGradient();
    Vector3 gradient = fullyFunctionalGradient.at(0);

    EXPECT_DOUBLE_EQ(gradient[0], 0.0);
    EXPECT_DOUBLE_EQ(gradient[1], 0.0);
    EXPECT_DOUBLE_EQ(gradient[2], 0.0);

    // Test CV constrained in x dimension
    std::vector<Vector3> constrainedInXGradient = constrainedInX->GetGradient();
    gradient = constrainedInXGradient.at(0);

    EXPECT_DOUBLE_EQ(gradient[0], 0.0);
    EXPECT_NEAR(gradient[1], -0.2922185356, eps);
    EXPECT_NEAR(gradient[2], -0.9563515711, eps);

    // Test fully constrained CV
    std::vector<Vector3> fullyConstrainedGradient = fullyConstrained->GetGradient();
    gradient = fullyConstrainedGradient.at(0);

    EXPECT_DOUBLE_EQ(gradient[0], 0.0);
    EXPECT_DOUBLE_EQ(gradient[1], 0.0);
    EXPECT_DOUBLE_EQ(gradient[2], 0.0);
}

TEST_F(ParticlePositionCVTest, BoundariesTest) {
    std::array<double, 2> bounds = fullyFunctional->GetBoundaries();

    EXPECT_DOUBLE_EQ(bounds[0], 0.0);
    EXPECT_DOUBLE_EQ(bounds[1], 0.0);
}

TEST_F(ParticlePositionCVTest, DifferenceTest) {
    // Test values
    double value1 = -15.3;
    double value2 = 0.0;
    double value3 = 4.2;
    double value4 = 1928437.3;

    // Initialize CVs
    fullyFunctional->Initialize(*snapshot1);
    constrainedInX->Initialize(*snapshot1);
    fullyConstrained->Initialize(*snapshot1);

    // Evaluate CVs
    fullyFunctional->Evaluate(*snapshot1);
    constrainedInX->Evaluate(*snapshot1);
    fullyConstrained->Evaluate(*snapshot1);

    // Test if CV values are as expected
    EXPECT_NEAR(fullyFunctional->GetDifference(value1), 17.8079872408, eps);
    EXPECT_NEAR(fullyFunctional->GetDifference(value2), 2.5079872408, eps);
    EXPECT_NEAR(fullyFunctional->GetDifference(value3), -1.6920127592, eps);
    EXPECT_NEAR(fullyFunctional->GetDifference(value4), -1928434.792, 0.001);

    EXPECT_NEAR(constrainedInX->GetDifference(value1), 19.0013511047, eps);
    EXPECT_NEAR(constrainedInX->GetDifference(value2), 3.7013511047, eps);
    EXPECT_NEAR(constrainedInX->GetDifference(value3), -0.4986488953, eps);
    EXPECT_NEAR(constrainedInX->GetDifference(value4), -1928433.599, 0.001);

    EXPECT_DOUBLE_EQ(fullyConstrained->GetDifference(value1), 15.3);
    EXPECT_DOUBLE_EQ(fullyConstrained->GetDifference(value2), 0.0);
    EXPECT_DOUBLE_EQ(fullyConstrained->GetDifference(value3), -4.2);
    EXPECT_DOUBLE_EQ(fullyConstrained->GetDifference(value4), -1928437.3);
}

TEST_F(ParticlePositionCVTest, SerializeTest) {
    // Implement Serialize Function first
}

int main(int argc, char *argv[])
{
    ::testing::InitGoogleTest(&argc, argv);
    boost::mpi::environment env(argc,argv);
    int ret = RUN_ALL_TESTS();
    return ret;
}

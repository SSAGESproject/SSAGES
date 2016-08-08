#include "../src/Utility/NearestNeighbor.h"
#include "gtest/gtest.h"

#define _USE_MATH_DEFINES
#include <cmath>

using namespace SSAGES;

// Test up to accuracy of 10^-10
const double eps = 1e-10;

// 90 degree angle
const double piOver2 = 0.5*M_PI;

class NearestNeighborTest : public ::testing::Test {

protected:
    virtual void SetUp() {
        cubic = {1.0, 1.0, 1.0, piOver2, piOver2, piOver2};
        orthorombic = {3.2, 1.4, 1.1, piOver2, piOver2, piOver2};
        monoclinic = {3.2, 1.4, 1.1, piOver2, 1.05, piOver2};
        triclinic = {3.2, 1.4, 1.1, 1.2, 1.05, 0.5};

        inBox1 = {0.2, 0.3, 0.9};
        inBox2 = {0.9, 0.7, 0.1};
        outOfBox1 = {4.2, -3.7, -0.1};
    }

    // Unit cells
    std::array<double, 6> cubic;
    std::array<double, 6> orthorombic;
    std::array<double, 6> monoclinic;
    std::array<double, 6> triclinic;

    // Vectors
    Vector3 inBox1;
    Vector3 inBox2;
    Vector3 outOfBox1;
};

TEST_F(NearestNeighborTest, samePositionTest)
{
    Vector3 vecCubic = NearestNeighbor(cubic, inBox1, inBox1);
    Vector3 vecOrthorombic = NearestNeighbor(orthorombic, inBox2, inBox2);
    Vector3 vecMonoclinic = NearestNeighbor(monoclinic, outOfBox1, outOfBox1);
    Vector3 vecTriclinic = NearestNeighbor(triclinic, inBox1, inBox1);

    // Test result, all Vectors should be zero
    // Cubic box
    EXPECT_NEAR(vecCubic[0], 0.0, eps);
    EXPECT_NEAR(vecCubic[1], 0.0, eps);
    EXPECT_NEAR(vecCubic[2], 0.0, eps);

    // Orthorombic box
    EXPECT_NEAR(vecOrthorombic[0], 0.0, eps);
    EXPECT_NEAR(vecOrthorombic[1], 0.0, eps);
    EXPECT_NEAR(vecOrthorombic[2], 0.0, eps);

    // Monoclinic box
    EXPECT_NEAR(vecMonoclinic[0], 0.0, eps);
    EXPECT_NEAR(vecMonoclinic[1], 0.0, eps);
    EXPECT_NEAR(vecMonoclinic[2], 0.0, eps);

    // Tricinlic box
    EXPECT_NEAR(vecTriclinic[0], 0.0, eps);
    EXPECT_NEAR(vecTriclinic[1], 0.0, eps);
    EXPECT_NEAR(vecTriclinic[2], 0.0, eps);
}

TEST_F(NearestNeighborTest, samePositionInMirrorBoxTest)
{
    Vector3 result = NearestNeighbor(cubic, inBox1, outOfBox1);

    EXPECT_NEAR(result[0], 0.0, eps);
    EXPECT_NEAR(result[1], 0.0, eps);
    EXPECT_NEAR(result[2], 0.0, eps);
}

TEST_F(NearestNeighborTest, cubicBoxTest)
{
    Vector3 result1 = NearestNeighbor(cubic, inBox1, inBox2);
    Vector3 result2 = NearestNeighbor(cubic, inBox2, outOfBox1);
    // inBox1 <-> outOfBox1 covered in samePositionInMirrorBoxTest

    EXPECT_NEAR(result1[0], 0.3, eps);
    EXPECT_NEAR(result1[1], -0.4, eps);
    EXPECT_NEAR(result1[2], -0.2, eps);

    EXPECT_NEAR(result2[0], -0.3, eps);
    EXPECT_NEAR(result2[1], 0.4, eps);
    EXPECT_NEAR(result2[2], 0.2, eps);
}

TEST_F(NearestNeighborTest, orthorombicBoxTest)
{
    Vector3 result1 = NearestNeighbor(orthorombic, inBox1, inBox2);
    Vector3 result2 = NearestNeighbor(orthorombic, inBox1, outOfBox1);
    Vector3 result3 = NearestNeighbor(orthorombic, inBox2, outOfBox1);

    EXPECT_NEAR(result1[0], -0.7, eps);
    EXPECT_NEAR(result1[1], -0.4, eps);
    EXPECT_NEAR(result1[2], -0.3, eps);

    EXPECT_NEAR(result2[0], -0.8, eps);
    EXPECT_NEAR(result2[1], -0.2, eps);
    EXPECT_NEAR(result2[2], -0.1, eps);

    EXPECT_NEAR(result3[0], -0.1, eps);
    EXPECT_NEAR(result3[1], 0.2, eps);
    EXPECT_NEAR(result3[2], 0.2, eps);
}

TEST_F(NearestNeighborTest, monoclinicBoxTest)
{
    Vector3 result1 = NearestNeighbor(monoclinic, inBox1, inBox2);
    Vector3 result2 = NearestNeighbor(monoclinic, inBox1, outOfBox1);
    Vector3 result3 = NearestNeighbor(monoclinic, inBox2, outOfBox1);

    EXPECT_NEAR(result1[0], -1.2473281527, eps);
    EXPECT_NEAR(result1[1], -0.4, eps);
    EXPECT_NEAR(result1[2], -0.1541655482, eps);

    EXPECT_NEAR(result2[0], -1.3473281527, eps);
    EXPECT_NEAR(result2[1], -0.2, eps);
    EXPECT_NEAR(result2[2], 0.0458344518, eps);

    EXPECT_NEAR(result3[0], -0.1, eps);
    EXPECT_NEAR(result3[1], 0.2, eps);
    EXPECT_NEAR(result3[2], 0.2, eps);
}

TEST_F(NearestNeighborTest, triclinicBoxTest)
{
    Vector3 result1 = NearestNeighbor(triclinic, inBox1, inBox2);
    Vector3 result2 = NearestNeighbor(triclinic, inBox1, outOfBox1);
    Vector3 result3 = NearestNeighbor(triclinic, inBox2, outOfBox1);

    EXPECT_NEAR(result1[0], -1.2473281527, eps);
    EXPECT_NEAR(result1[1], -0.2295207370, eps);
    EXPECT_NEAR(result1[2], -0.1388123956, eps);

    EXPECT_NEAR(result2[0], 0.8809783274, eps);
    EXPECT_NEAR(result2[1], 0.1433047388, eps);
    EXPECT_NEAR(result2[2], 0.0611876044, eps);

    EXPECT_NEAR(result3[0], 0.8996908935, eps);
    EXPECT_NEAR(result3[1], -0.2983702783, eps);
    EXPECT_NEAR(result3[2], 0.2, eps);
}

int main(int argc, char *argv[])
{
    ::testing::InitGoogleTest(&argc, argv);
    int ret = RUN_ALL_TESTS();
    return ret;
}

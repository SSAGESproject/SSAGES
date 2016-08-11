#include "../src/Utility/UnwrapCoordinates.h"
#include "gtest/gtest.h"

#define _USE_MATH_DEFINES
#include <cmath>

using namespace SSAGES;

// Test up to accuracy of 10^-10
const double eps = 1e-10;

// 90 degree angle
const double piOver2 = 0.5*M_PI;

class UnwrapCoordinatesTest : public ::testing::Test {

protected:
    virtual void SetUp() {
        cubic = {1.0, 1.0, 1.0, piOver2, piOver2, piOver2};
        orthorombic = {3.2, 1.4, 1.1, piOver2, piOver2, piOver2};
        monoclinic = {3.2, 1.4, 1.1, piOver2, 1.05, piOver2};
        triclinic = {3.2, 1.4, 1.1, 1.2, 1.05, 0.5};

        inBox1 = {0.2, 0.3, 0.9};
        inBox2 = {0.9, 0.7, 0.1};
        outOfBox1 = {4.2, -3.7, -0.1};

        inBox1Mirrors = {0, 0, 0};
        inBox2Mirrors = {1, 1, -1};
        outOfBox1Mirrors = {1, -2, 0};
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

    // Corresponding Mirror images
    Vector3 inBox1Mirrors;
    Vector3 inBox2Mirrors;
    Vector3 outOfBox1Mirrors;
};

TEST_F(UnwrapCoordinatesTest, cubicBoxTest)
{
    // Test cubic box (1.0, 1.0, 1.0), (piOver2, piOver2, piOver2)
    Vector3 result1 = UnwrapCoordinates(cubic, inBox1, inBox1Mirrors);
    Vector3 result2 = UnwrapCoordinates(cubic, inBox2, inBox2Mirrors);
    Vector3 result3 = UnwrapCoordinates(cubic, outOfBox1, outOfBox1Mirrors);

    // Test inBox1 (0.2, 0.3, 0.9), mirror (0, 0, 0)
    EXPECT_NEAR(result1[0], 0.2, eps);
    EXPECT_NEAR(result1[1], 0.3, eps);
    EXPECT_NEAR(result1[2], 0.9, eps);

    // Test inBox2 (0.9, 0.7, 0.1), mirror (1, 1, -1)
    EXPECT_NEAR(result2[0], 1.9, eps);
    EXPECT_NEAR(result2[1], 1.7, eps);
    EXPECT_NEAR(result2[2], -0.9, eps);

    // Test outOfBox1 (4.2, -3.7, -0.1), mirror (1, -2, 0)
    EXPECT_NEAR(result3[0], 5.2, eps);
    EXPECT_NEAR(result3[1], -5.7, eps);
    EXPECT_NEAR(result3[2], -0.1, eps);
}

TEST_F(UnwrapCoordinatesTest, orthorombicBoxTest)
{
    // Test orthorombic box (3.2, 1.4, 1.1), (piOver2, piOver2, piOver2)
    Vector3 result1 = UnwrapCoordinates(orthorombic, inBox1, inBox1Mirrors);
    Vector3 result2 = UnwrapCoordinates(orthorombic, inBox2, inBox2Mirrors);
    Vector3 result3 = UnwrapCoordinates(orthorombic, outOfBox1, outOfBox1Mirrors);

    // Test inBox1 (0.2, 0.3, 0.9), mirror (0, 0, 0)
    EXPECT_NEAR(result1[0], 0.2, eps);
    EXPECT_NEAR(result1[1], 0.3, eps);
    EXPECT_NEAR(result1[2], 0.9, eps);

    // Test inBox2 (0.9, 0.7, 0.1), mirror (1, 1, -1)
    EXPECT_NEAR(result2[0], 4.1, eps);
    EXPECT_NEAR(result2[1], 2.1, eps);
    EXPECT_NEAR(result2[2], -1.0, eps);

    // Test outOfBox1 (4.2, -3.7, -0.1), mirror (1, -2, 0)
    EXPECT_NEAR(result3[0], 7.4, eps);
    EXPECT_NEAR(result3[1], -6.5, eps);
    EXPECT_NEAR(result3[2], -0.1, eps);
}

TEST_F(UnwrapCoordinatesTest, monoclinicBoxTest)
{
    // Test monoclinic box (3.2, 1.4, 1.1), (piOver2, 1.05, piOver2)
    Vector3 result1 = UnwrapCoordinates(monoclinic, inBox1, inBox1Mirrors);
    Vector3 result2 = UnwrapCoordinates(monoclinic, inBox2, inBox2Mirrors);
    Vector3 result3 = UnwrapCoordinates(monoclinic, outOfBox1, outOfBox1Mirrors);

    // Test inBox1 (0.2, 0.3, 0.9), mirror (0, 0, 0)
    EXPECT_NEAR(result1[0], 0.2, eps);
    EXPECT_NEAR(result1[1], 0.3, eps);
    EXPECT_NEAR(result1[2], 0.9, eps);

    // Test inBox2 (0.9, 0.7, 0.1), mirror (1, 1, -1)
    EXPECT_NEAR(result2[0], 3.5526718473, eps);
    EXPECT_NEAR(result2[1], 2.1, eps);
    EXPECT_NEAR(result2[2], -0.8541655482, eps);

    // Test outOfBox1 (4.2, -3.7, -0.1), mirror (1, -2, 0)
    EXPECT_NEAR(result3[0], 7.4, eps);
    EXPECT_NEAR(result3[1], -6.5, eps);
    EXPECT_NEAR(result3[2], -0.1, eps);
}

TEST_F(UnwrapCoordinatesTest, triclinicBoxTest)
{
    // Test triclinic box (3.2, 1.4, 1.1) (1.2, 1.05, 0.5)
    Vector3 result1 = UnwrapCoordinates(triclinic, inBox1, inBox1Mirrors);
    Vector3 result2 = UnwrapCoordinates(triclinic, inBox2, inBox2Mirrors);
    Vector3 result3 = UnwrapCoordinates(triclinic, outOfBox1, outOfBox1Mirrors);

    // Test inBox1 (0.2, 0.3, 0.9), mirror (0, 0, 0)
    EXPECT_NEAR(result1[0], 0.2, eps);
    EXPECT_NEAR(result1[1], 0.3, eps);
    EXPECT_NEAR(result1[2], 0.9, eps);

    // Test inBox2 (0.9, 0.7, 0.1), mirror (1, 1, -1)
    EXPECT_NEAR(result2[0], 4.7812874340, eps);
    EXPECT_NEAR(result2[1], 1.5416750171, eps);
    EXPECT_NEAR(result2[2], -0.8388123956, eps);

    // Test outOfBox1 (4.2, -3.7, -0.1), mirror (1, -2, 0)
    EXPECT_NEAR(result3[0], 4.9427688267, eps);
    EXPECT_NEAR(result3[1], -5.0423915081, eps);
    EXPECT_NEAR(result3[2], -0.1, eps);
}

int main(int argc, char *argv[])
{
    ::testing::InitGoogleTest(&argc, argv);
    int ret = RUN_ALL_TESTS();
    return ret;
}
#include "../src/Utility/UnitCellConversion.h"
#include "gtest/gtest.h"

#define _USE_MATH_DEFINES
#include <cmath>
#include <stdexcept>

using namespace SSAGES;

// Test up to accuracy of 10^-10
const double eps = 1e-10;

// 90 degree angle
const double piOver2 = 0.5*M_PI;

class UnitCellConversionTest : public ::testing::Test {

protected:
    virtual void SetUp() {
        cubic = {1.0, 1.0, 1.0, piOver2, piOver2, piOver2};
        orthorombic = {3.2, 1.4, 1.1, piOver2, piOver2, piOver2};
        monoclinic = {3.2, 1.4, 1.1, piOver2, 1.05, piOver2};
        triclinic = {3.2, 1.4, 1.1, 1.2, 1.05, 0.5};
    }

    std::array<double, 6> cubic;
    std::array<double, 6> orthorombic;
    std::array<double, 6> monoclinic;
    std::array<double, 6> triclinic;
};

TEST_F(UnitCellConversionTest, IllegalArgumentsTest)
{
    std::array<double, 6> negativeX = {-2.3, 1.0, 1.2, piOver2, piOver2, piOver2};
    std::array<double, 6> negativeY = {1.2, -1.1, 0.8, piOver2, piOver2, piOver2};
    std::array<double, 6> negativeZ = {0.3, 2.1, -1.1, piOver2, piOver2, piOver2};

    std::array<double, 6> negativeAlpha = {1.0, 1.0, 1.0, -0.34, piOver2, piOver2};
    std::array<double, 6> alphaLargerThanPi = {1.0, 1.0, 1.0, 3.42, piOver2, piOver2};

    std::array<double, 6> negativeBeta = {1.0, 1.0, 1.0, piOver2, -0.1, piOver2};
    std::array<double, 6> betaLargerThanPi = {1.0, 1.0, 1.0, piOver2, 4.1, piOver2};

    std::array<double, 6> negativeGamma = {1.0, 1.0, 1.0, piOver2, piOver2, -0.3};
    std::array<double, 6> gammaLargerThanPi = {1.0, 1.0, 1.0, piOver2, piOver2, 3.15};

    // Test uc_c2f
    EXPECT_THROW(uc_c2f(negativeX), std::domain_error);
    EXPECT_THROW(uc_c2f(negativeY), std::domain_error);
    EXPECT_THROW(uc_c2f(negativeZ), std::domain_error);

    EXPECT_THROW(uc_c2f(negativeAlpha), std::domain_error);
    EXPECT_THROW(uc_c2f(alphaLargerThanPi), std::domain_error);

    EXPECT_THROW(uc_c2f(negativeBeta), std::domain_error);
    EXPECT_THROW(uc_c2f(betaLargerThanPi), std::domain_error);

    EXPECT_THROW(uc_c2f(negativeGamma), std::domain_error);
    EXPECT_THROW(uc_c2f(gammaLargerThanPi), std::domain_error);

    // Test uc_f2c
    EXPECT_THROW(uc_f2c(negativeX), std::domain_error);
    EXPECT_THROW(uc_f2c(negativeY), std::domain_error);
    EXPECT_THROW(uc_f2c(negativeZ), std::domain_error);

    EXPECT_THROW(uc_f2c(negativeAlpha), std::domain_error);
    EXPECT_THROW(uc_f2c(alphaLargerThanPi), std::domain_error);

    EXPECT_THROW(uc_f2c(negativeBeta), std::domain_error);
    EXPECT_THROW(uc_f2c(betaLargerThanPi), std::domain_error);

    EXPECT_THROW(uc_f2c(negativeGamma), std::domain_error);
    EXPECT_THROW(uc_f2c(gammaLargerThanPi), std::domain_error);
}

TEST_F(UnitCellConversionTest, cartesian2FractionalTest)
{
    // Test cubic cell (1.0, 1.0, 1.0, piOver2, piOver2, piOver2)
    std::array<std::array<double,3>,3> cubic_matrix = uc_c2f(cubic);

    EXPECT_NEAR(cubic_matrix[0][0], 1.0, eps);
    EXPECT_NEAR(cubic_matrix[0][1], 0.0, eps);
    EXPECT_NEAR(cubic_matrix[0][2], 0.0, eps);

    EXPECT_NEAR(cubic_matrix[1][0], 0.0, eps);
    EXPECT_NEAR(cubic_matrix[1][1], 1.0, eps);
    EXPECT_NEAR(cubic_matrix[1][2], 0.0, eps);

    EXPECT_NEAR(cubic_matrix[2][0], 0.0, eps);
    EXPECT_NEAR(cubic_matrix[2][1], 0.0, eps);
    EXPECT_NEAR(cubic_matrix[2][2], 1.0, eps);

    // Test orthorombic cell (3.2, 1.4, 1.1, piOver2, piOver2, piOver2)
    std::array<std::array<double,3>,3> orthorombic_matrix = uc_c2f(orthorombic);

    EXPECT_NEAR(orthorombic_matrix[0][0], 0.3125, eps);
    EXPECT_NEAR(orthorombic_matrix[0][1], 0.0, eps);
    EXPECT_NEAR(orthorombic_matrix[0][2], 0.0, eps);

    EXPECT_NEAR(orthorombic_matrix[1][0], 0.0, eps);
    EXPECT_NEAR(orthorombic_matrix[1][1], 0.7142857143, eps);
    EXPECT_NEAR(orthorombic_matrix[1][2], 0.0, eps);

    EXPECT_NEAR(orthorombic_matrix[2][0], 0.0, eps);
    EXPECT_NEAR(orthorombic_matrix[2][1], 0.0, eps);
    EXPECT_NEAR(orthorombic_matrix[2][2], 0.9090909091, eps);

    // Test monoclinic cell (3.2, 1.4, 1.1, piOver2, 1.05, piOver2)
    std::array<std::array<double,3>,3> monoclinic_matrix = uc_c2f(monoclinic);

    EXPECT_NEAR(monoclinic_matrix[0][0], 0.3125, eps);
    EXPECT_NEAR(monoclinic_matrix[0][1], 0.0, eps);
    EXPECT_NEAR(monoclinic_matrix[0][2], -0.1792561553, eps);

    EXPECT_NEAR(monoclinic_matrix[1][0], 0.0, eps);
    EXPECT_NEAR(monoclinic_matrix[1][1], 0.7142857143, eps);
    EXPECT_NEAR(monoclinic_matrix[1][2], 0.0, eps);

    EXPECT_NEAR(monoclinic_matrix[2][0], 0.0, eps);
    EXPECT_NEAR(monoclinic_matrix[2][1], 0.0, eps);
    EXPECT_NEAR(monoclinic_matrix[2][2], 1.0480361630, eps);

    // Test triclinic cell (3.2, 1.4, 1.1, 1.2, 1.05, 0.5)
    std::array<std::array<double,3>,3> triclinic_matrix = uc_c2f(triclinic);

    EXPECT_NEAR(triclinic_matrix[0][0], 0.3125, eps);
    EXPECT_NEAR(triclinic_matrix[0][1], -0.5720274130, eps);
    EXPECT_NEAR(triclinic_matrix[0][2], -0.2860623281, eps);

    EXPECT_NEAR(triclinic_matrix[1][0], 0.0, eps);
    EXPECT_NEAR(triclinic_matrix[1][1], 1.4898783164, eps);
    EXPECT_NEAR(triclinic_matrix[1][2], 0.2705475115, eps);

    EXPECT_NEAR(triclinic_matrix[2][0], 0.0, eps);
    EXPECT_NEAR(triclinic_matrix[2][1], 0.0, eps);
    EXPECT_NEAR(triclinic_matrix[2][2], 1.0651755395, eps);
}

TEST_F(UnitCellConversionTest, fractional2CartesianTest)
{
    // Test cubic cell (1.0, 1.0, 1.0, piOver2, piOver2, piOver2)
    std::array<std::array<double,3>,3> cubic_matrix = uc_f2c(cubic);

    EXPECT_NEAR(cubic_matrix[0][0], 1.0, eps);
    EXPECT_NEAR(cubic_matrix[0][1], 0.0, eps);
    EXPECT_NEAR(cubic_matrix[0][2], 0.0, eps);

    EXPECT_NEAR(cubic_matrix[1][0], 0.0, eps);
    EXPECT_NEAR(cubic_matrix[1][1], 1.0, eps);
    EXPECT_NEAR(cubic_matrix[1][2], 0.0, eps);

    EXPECT_NEAR(cubic_matrix[2][0], 0.0, eps);
    EXPECT_NEAR(cubic_matrix[2][1], 0.0, eps);
    EXPECT_NEAR(cubic_matrix[2][2], 1.0, eps);

    // Test orthorombic cell (3.2, 1.4, 1.1, piOver2, piOver2, piOver2)
    std::array<std::array<double,3>,3> orthorombic_matrix = uc_f2c(orthorombic);

    EXPECT_NEAR(orthorombic_matrix[0][0], 3.2, eps);
    EXPECT_NEAR(orthorombic_matrix[0][1], 0.0, eps);
    EXPECT_NEAR(orthorombic_matrix[0][2], 0.0, eps);

    EXPECT_NEAR(orthorombic_matrix[1][0], 0.0, eps);
    EXPECT_NEAR(orthorombic_matrix[1][1], 1.4, eps);
    EXPECT_NEAR(orthorombic_matrix[1][2], 0.0, eps);

    EXPECT_NEAR(orthorombic_matrix[2][0], 0.0, eps);
    EXPECT_NEAR(orthorombic_matrix[2][1], 0.0, eps);
    EXPECT_NEAR(orthorombic_matrix[2][2], 1.1, eps);

    // Test monoclinic cell (3.2, 1.4, 1.1, piOver2, 1.05, piOver2)
    std::array<std::array<double,3>,3> monoclinic_matrix = uc_f2c(monoclinic);

    EXPECT_NEAR(monoclinic_matrix[0][0], 3.2, eps);
    EXPECT_NEAR(monoclinic_matrix[0][1], 0.0, eps);
    EXPECT_NEAR(monoclinic_matrix[0][2], 0.5473281527, eps);

    EXPECT_NEAR(monoclinic_matrix[1][0], 0.0, eps);
    EXPECT_NEAR(monoclinic_matrix[1][1], 1.4, eps);
    EXPECT_NEAR(monoclinic_matrix[1][2], 0.0, eps);

    EXPECT_NEAR(monoclinic_matrix[2][0], 0.0, eps);
    EXPECT_NEAR(monoclinic_matrix[2][1], 0.0, eps);
    EXPECT_NEAR(monoclinic_matrix[2][2], 0.9541655482, eps);

    // Test triclinic cell (3.2, 1.4, 1.1, 1.2, 1.05, 0.5)
    std::array<std::array<double,3>,3> triclinic_matrix = uc_f2c(triclinic);

    EXPECT_NEAR(triclinic_matrix[0][0], 3.2, eps);
    EXPECT_NEAR(triclinic_matrix[0][1], 1.2286155866, eps);
    EXPECT_NEAR(triclinic_matrix[0][2], 0.5473281527, eps);

    EXPECT_NEAR(triclinic_matrix[1][0], 0.0, eps);
    EXPECT_NEAR(triclinic_matrix[1][1], 0.6711957540, eps);
    EXPECT_NEAR(triclinic_matrix[1][2], -0.1704792630, eps);

    EXPECT_NEAR(triclinic_matrix[2][0], 0.0, eps);
    EXPECT_NEAR(triclinic_matrix[2][1], 0.0, eps);
    EXPECT_NEAR(triclinic_matrix[2][2], 0.9388123956, eps);
}

int main(int argc, char *argv[])
{
    ::testing::InitGoogleTest(&argc, argv);
    int ret = RUN_ALL_TESTS();
    return ret;
}

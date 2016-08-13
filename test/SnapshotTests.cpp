#include "../src/Snapshot.h"
#include "gtest/gtest.h"

using namespace SSAGES;

// Test up to accuracy of 10^-10
const double eps = 1e-10;

class SnapshotTest : public ::testing::Test 
{
protected:
    boost::mpi::communicator comm;
	
    // Snapshots.
	std::shared_ptr<Snapshot> cubic, orthorhombic, monoclinic, triclinic; 

	// Vectors.
	Vector3 inBox1, inBox2, outOfBox1;

	virtual void SetUp()
	{
		Matrix3 H; 
		cubic = std::make_shared<Snapshot>(comm, 0);
		H << 1.0, 0.0, 0.0,
		     0.0, 1.0, 0.0,
		     0.0, 0.0, 1.0; 
		cubic->SetHMatrix(H);

		orthorhombic = std::make_shared<Snapshot>(comm, 0);
		H << 3.2, 0.0, 0.0,
		     0.0, 1.4, 0.0,
		     0.0, 0.0, 1.1; 
		orthorhombic->SetHMatrix(H);

		monoclinic = std::make_shared<Snapshot>(comm, 0);
		H << 3.2000, 0.0000, 0.547328152680900,
		          0, 1.4000,            0.0000,
		          0,      0, 0.954165548153419;
		monoclinic->SetHMatrix(H);

		triclinic = std::make_shared<Snapshot>(comm, 0);
		H << 3.2000, 1.228615586646522,  0.547328152680900,
		          0, 0.671195754045884, -0.170479263032258,
		          0,                 0,  0.938812395614210;
		triclinic->SetHMatrix(H);

		inBox1 = {0.2, 0.3, 0.9};
        inBox2 = {0.9, 0.7, 0.1};
        outOfBox1 = {4.2, -3.7, -0.1};
	}
};

TEST_F(SnapshotTest, SamePositionTest)
{
	Vector3 vCubic = inBox1 - inBox1;
	cubic->ApplyMinimumImage(&vCubic);

	Vector3 vOrtho = inBox2 - inBox2;
	orthorhombic->ApplyMinimumImage(&vOrtho);

	Vector3 vMonoc = outOfBox1 - outOfBox1;
	monoclinic->ApplyMinimumImage(&vMonoc);

	Vector3 vTricl = inBox1 - inBox1;
	triclinic->ApplyMinimumImage(&vTricl);

    // Test result, all Vectors should be zero
    // Cubic box
    EXPECT_NEAR(vCubic[0], 0.0, eps);
    EXPECT_NEAR(vCubic[1], 0.0, eps);
    EXPECT_NEAR(vCubic[2], 0.0, eps);

    // Orthorombic box
    EXPECT_NEAR(vOrtho[0], 0.0, eps);
    EXPECT_NEAR(vOrtho[1], 0.0, eps);
    EXPECT_NEAR(vOrtho[2], 0.0, eps);

    // Monoclinic box
    EXPECT_NEAR(vMonoc[0], 0.0, eps);
    EXPECT_NEAR(vMonoc[1], 0.0, eps);
    EXPECT_NEAR(vMonoc[2], 0.0, eps);

    // Tricinlic box
    EXPECT_NEAR(vTricl[0], 0.0, eps);
    EXPECT_NEAR(vTricl[1], 0.0, eps);
    EXPECT_NEAR(vTricl[2], 0.0, eps);
}

TEST_F(SnapshotTest, SamePositionInMirrorBoxTest)
{
	Vector3 dx = inBox1 - outOfBox1;
	cubic->ApplyMinimumImage(&dx);
	
	EXPECT_NEAR(dx[0], 0.0, eps);
    EXPECT_NEAR(dx[1], 0.0, eps);
    EXPECT_NEAR(dx[2], 0.0, eps);
}

TEST_F(SnapshotTest, CubicBoxTest)
{
	Vector3 dx1 = inBox1 - inBox2; 
	Vector3 dx2 = inBox2 - outOfBox1;
    // inBox1 <-> outOfBox1 covered in samePositionInMirrorBoxTest
    cubic->ApplyMinimumImage(&dx1);
    cubic->ApplyMinimumImage(&dx2);

	EXPECT_NEAR(dx1[0], 0.3, eps);
	EXPECT_NEAR(dx1[1], -0.4, eps);
	EXPECT_NEAR(dx1[2], -0.2, eps);

	EXPECT_NEAR(dx2[0], -0.3, eps);
	EXPECT_NEAR(dx2[1], 0.4, eps);
	EXPECT_NEAR(dx2[2], 0.2, eps);
}

TEST_F(SnapshotTest, OrthorhombicBoxTest)
{
	Vector3 dx1 = inBox1 - inBox2;
	Vector3 dx2 = inBox1 - outOfBox1;
	Vector3 dx3 = inBox2 - outOfBox1;
	orthorhombic->ApplyMinimumImage(&dx1);
	orthorhombic->ApplyMinimumImage(&dx2);
	orthorhombic->ApplyMinimumImage(&dx3);

	EXPECT_NEAR(dx1[0], -0.7, eps);
	EXPECT_NEAR(dx1[1], -0.4, eps);
	EXPECT_NEAR(dx1[2], -0.3, eps);

	EXPECT_NEAR(dx2[0], -0.8, eps);
	EXPECT_NEAR(dx2[1], -0.2, eps);
	EXPECT_NEAR(dx2[2], -0.1, eps);

	EXPECT_NEAR(dx3[0], -0.1, eps);
	EXPECT_NEAR(dx3[1], 0.2, eps);
	EXPECT_NEAR(dx3[2], 0.2, eps);
}

TEST_F(SnapshotTest, MonoclinicBoxTest)
{
	Vector3 dx1 = inBox1 - inBox2;
	Vector3 dx2 = inBox1 - outOfBox1;
	Vector3 dx3 = inBox2 - outOfBox1;
	monoclinic->ApplyMinimumImage(&dx1);
	monoclinic->ApplyMinimumImage(&dx2);
	monoclinic->ApplyMinimumImage(&dx3);
    
    EXPECT_NEAR(dx1[0], -1.2473281527, eps);
    EXPECT_NEAR(dx1[1], -0.4, eps);
    EXPECT_NEAR(dx1[2], -0.1541655482, eps);

    EXPECT_NEAR(dx2[0], -1.3473281527, eps);
    EXPECT_NEAR(dx2[1], -0.2, eps);
    EXPECT_NEAR(dx2[2], 0.0458344518, eps);

    EXPECT_NEAR(dx3[0], -0.1, eps);
    EXPECT_NEAR(dx3[1], 0.2, eps);
    EXPECT_NEAR(dx3[2], 0.2, eps);
}

TEST_F(SnapshotTest, TriclinicBoxTest)
{
	Vector3 dx1 = inBox1 - inBox2;
	Vector3 dx2 = inBox1 - outOfBox1;
	Vector3 dx3 = inBox2 - outOfBox1;
	triclinic->ApplyMinimumImage(&dx1);
	triclinic->ApplyMinimumImage(&dx2);
	triclinic->ApplyMinimumImage(&dx3);
    
	EXPECT_NEAR(dx1[0], -1.2473281527, eps);
	EXPECT_NEAR(dx1[1], -0.2295207370, eps);
	EXPECT_NEAR(dx1[2], -0.1388123956, eps);

	EXPECT_NEAR(dx2[0], 0.8809783274, eps);
	EXPECT_NEAR(dx2[1], 0.1433047388, eps);
	EXPECT_NEAR(dx2[2], 0.0611876044, eps);

	EXPECT_NEAR(dx3[0], 0.8996908935, eps);
	EXPECT_NEAR(dx3[1], -0.2983702783, eps);
	EXPECT_NEAR(dx3[2], 0.2, eps);
}

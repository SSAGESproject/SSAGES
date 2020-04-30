#include "Tests.h"

using namespace SSAGES;

class SnapshotTests : public ::testing::Test
{
protected:
    mxx::comm comm;
	
    // Snapshots.
	std::shared_ptr<Snapshot> cubic, orthorhombic, monoclinic, triclinic; 

	// Vectors.
	Vector3 inBox1, inBox2, outOfBox1;

	// Origins.
	Vector3 origin1, origin2;

	// Image flags (mirrors).
	Integer3 inBox1Mirrors, inBox2Mirrors, outOfBox1Mirrors;

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

        inBox1Mirrors = {0, 0, 0};
        inBox2Mirrors = {1, 1, -1};
        outOfBox1Mirrors = {1, -2, 0};

        origin1 = {0, 0, 0};
        origin2 = {-0.5, -0.3, -0.1};
	}
};

TEST_F(SnapshotTests, SamePositionTest)
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

TEST_F(SnapshotTests, SamePositionInMirrorBoxTest)
{
	Vector3 dx = inBox1 - outOfBox1;
	cubic->ApplyMinimumImage(&dx);
	
	EXPECT_NEAR(dx[0], 0.0, eps);
    EXPECT_NEAR(dx[1], 0.0, eps);
    EXPECT_NEAR(dx[2], 0.0, eps);
}

TEST_F(SnapshotTests, CubicBoxTest)
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

TEST_F(SnapshotTests, OrthorhombicBoxTest)
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

TEST_F(SnapshotTests, MonoclinicBoxTest)
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

TEST_F(SnapshotTests, TriclinicBoxTest)
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

TEST_F(SnapshotTests, UnwrapCubicBoxTest)
{
	auto result1 = cubic->UnwrapVector(inBox1, inBox1Mirrors);
	auto result2 = cubic->UnwrapVector(inBox2, inBox2Mirrors);
	auto result3 = cubic->UnwrapVector(outOfBox1, outOfBox1Mirrors);
	
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

TEST_F(SnapshotTests, UnwrapOrthorhombicBoxTest)
{
	auto result1 = orthorhombic->UnwrapVector(inBox1, inBox1Mirrors);
	auto result2 = orthorhombic->UnwrapVector(inBox2, inBox2Mirrors);
	auto result3 = orthorhombic->UnwrapVector(outOfBox1, outOfBox1Mirrors);
	
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

TEST_F(SnapshotTests, UnwrapMonoclinicBoxTest)
{
	auto result1 = monoclinic->UnwrapVector(inBox1, inBox1Mirrors);
	auto result2 = monoclinic->UnwrapVector(inBox2, inBox2Mirrors);
	auto result3 = monoclinic->UnwrapVector(outOfBox1, outOfBox1Mirrors);
	
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

TEST_F(SnapshotTests, UnwrapTriclinicBoxTest)
{
	auto result1 = triclinic->UnwrapVector(inBox1, inBox1Mirrors);
	auto result2 = triclinic->UnwrapVector(inBox2, inBox2Mirrors);
	auto result3 = triclinic->UnwrapVector(outOfBox1, outOfBox1Mirrors);

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

TEST_F(SnapshotTests, OriginCubicBoxTest)
{
	cubic->SetOrigin(origin1);
	auto result1 = cubic->ScaleVector(inBox1);
	EXPECT_NEAR(result1[0], 0.2, eps);
	EXPECT_NEAR(result1[1], 0.3, eps);
	EXPECT_NEAR(result1[2], 0.9, eps);

	cubic->SetOrigin(origin2);
	auto result2 = cubic->ScaleVector(inBox1);
	EXPECT_NEAR(result2[0], 0.7, eps);
	EXPECT_NEAR(result2[1], 0.6, eps);
	EXPECT_NEAR(result2[2], 1.0, eps);
}

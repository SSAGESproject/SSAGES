#include "../src/Snapshot.h"
#include "gtest/gtest.h"

using namespace SSAGES;

// Test up to accuracy of 10^-10
constexpr double eps = 1e-10;

class COMTest : public ::testing::Test
{
protected:
	boost::mpi::communicator comm;

    // Snapshots.
	std::shared_ptr<Snapshot> cubic;

	// Positions. 
	std::vector<Vector3> cube3d;

	// Indices.
	Label indices;

	// Masses.
	std::vector<double> masses;

	virtual void SetUp()
	{
		Matrix3 H;
		cubic = std::make_shared<Snapshot>(comm, 0);
		H << 10.0,    0,    0,
		        0, 10.0,    0,
		        0,    0, 10.0;
		cubic->SetHMatrix(H);

		if(comm.size() == 2)
		{
			if(comm.rank() == 0)
			{
				cube3d = {
					{8,     8,     8},
					{8,     8,     2},
					{8,     2,     8},
					{8,     2,     2},
				};
			}
			else
			{
				cube3d = {
					{ 2,     8,     8},
					{ 2,     8,     2},
					{ 2,     2,     8},
					{ 2,     2,     2}
				};
			}
			indices = {0, 1, 2, 3};
			masses = {1.0, 1.0, 1.0, 1.0};
		}
		else
		{
			cube3d = {
				{8,    8,    8},
				{8,    8,     2},
				{8,     2,    8},
				{8,     2,     2},
				{ 2,    8,    8},
				{ 2,    8,     2},
				{ 2,     2,    8},
				{ 2,     2,     2}
			};
			indices = {0, 1, 2, 3, 4, 5, 6, 7};
			masses = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
		}

		auto& pos = cubic->GetPositions();
		for(auto& c : cube3d)
			pos.push_back(c);

		auto& m = cubic->GetMasses();
		m = masses;
	}
};

TEST_F(COMTest, ZeroCOMTest)
{
	// The cube is sitting at a periodic corner, 
	// but COM should be at the origin (or opposite corner).
	auto com = cubic->CenterOfMass(indices);
	EXPECT_NEAR(com[0], 10.0, eps);
	EXPECT_NEAR(com[1], 10.0, eps);
	EXPECT_NEAR(com[2], 10.0, eps);
}

int main(int argc, char *argv[])
{
    ::testing::InitGoogleTest(&argc, argv);
     boost::mpi::environment env(argc,argv);
    int ret = RUN_ALL_TESTS();

    return ret;
}
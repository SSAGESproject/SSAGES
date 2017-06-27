#include "../src/CVs/NearestNeighborsCV.h"
#include "../src/Snapshot.h"
#include "gtest/gtest.h"
#include <mxx/env.hpp>
#include <mxx/comm.hpp>

using namespace SSAGES;

// Test up to accuracy of 10^-10
constexpr double eps = 1e-10;

class NearestNeighborsCVTest : public ::testing::Test{
protected:
	virtual void SetUp()
	{
		cv1 = new NearestNeighborsCV({1, 2, 3, 4}, GaussianFunction(1.0, 1.0));
		cv2 = new NearestNeighborsCV({1, 2, 3, 4}, GaussianFunction(2.0, 1.0));

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
	NearestNeighborsCV *cv1, *cv2;
	mxx::comm comm;
};

TEST_F(NearestNeighborsCVTest, DefaultBehavior)
{
	cv1->Initialize(*snapshot);
	cv1->Evaluate(*snapshot);
	EXPECT_NEAR(cv1->GetValue(), 5.300240013212994, eps);

	auto& grad1 = cv1->GetGradient(); 
	// NOTE: Gradient checks are only valid for single processor
	// but value checks for 2. This needs to be adjusted later. 
	// The gradients ARE correct for 2 cores, but the global/local
	// indices don't match... to be fixed later.
	EXPECT_NEAR(grad1[0][0],  0.375706014812775, eps);
	EXPECT_NEAR(grad1[0][1], -0.320945252301533, eps);
	EXPECT_NEAR(grad1[0][2],  0.534999020276341, eps);

	EXPECT_NEAR(grad1[1][0], -0.134120305713599, eps);
	EXPECT_NEAR(grad1[1][1],  0.248756529603423, eps);
	EXPECT_NEAR(grad1[1][2], -0.184465705213535, eps);

	EXPECT_NEAR(grad1[2][0], -0.179139027312953, eps);
	EXPECT_NEAR(grad1[2][1], -0.089463524139856, eps);
	EXPECT_NEAR(grad1[2][2], -0.089675503173096, eps);

	EXPECT_NEAR(grad1[3][0], -0.062446681786223, eps);
	EXPECT_NEAR(grad1[3][1],  0.161652246837967, eps);
	EXPECT_NEAR(grad1[3][2], -0.260857811889710, eps);

	
	cv2->Initialize(*snapshot);
	cv2->Evaluate(*snapshot);

	EXPECT_NEAR(cv2->GetValue(), 4.579273659770989, eps);

	auto& grad2 = cv2->GetGradient();

	EXPECT_NEAR(grad2[0][0], -0.266194699985101, eps);
	EXPECT_NEAR(grad2[0][1],  0.108858442919108, eps);
	EXPECT_NEAR(grad2[0][2], -0.300430374114451, eps);
	
	EXPECT_NEAR(grad2[1][0], -0.298226633238643, eps);
	EXPECT_NEAR(grad2[1][1], -0.268582500391991, eps);
	EXPECT_NEAR(grad2[1][2],  0.366697981497343, eps);
	
	EXPECT_NEAR(grad2[2][0],  0.291627507261988, eps);
	EXPECT_NEAR(grad2[2][1],  0.537612156118958, eps);
	EXPECT_NEAR(grad2[2][2],  0.057280680999346, eps);
	
	EXPECT_NEAR(grad2[3][0],  0.272793825961756, eps);
	EXPECT_NEAR(grad2[3][1], -0.377888098646075, eps);
	EXPECT_NEAR(grad2[3][2], -0.123548288382238, eps);
}

int main(int argc, char *argv[])
{
	printf("Got to here.\n");
	::testing::InitGoogleTest(&argc, argv);
	mxx::env e(argc,argv);
	int ret = RUN_ALL_TESTS();
	return ret;
}
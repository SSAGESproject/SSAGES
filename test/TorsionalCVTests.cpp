#include "../src/CVs/TorsionalCV.h"
#include "../src/Snapshot.h"
#include "gtest/gtest.h"
#include <boost/mpi.hpp>

using namespace SSAGES;

TEST(TorsionalCV, DefaultBehavior)
{
	// Initialize atoms and CV.
	boost::mpi::communicator comm;

	auto snap = Snapshot(comm, 0);

	auto& pos = snap.GetPositions();
	pos.resize(6);
	auto& ids = snap.GetAtomIDs();
	ids.resize(6);

	for(unsigned int i =0; i <ids.size();i++)
		ids[i] = i+1;

	pos[0][0] = 0; 
	pos[0][1] = 0; 
	pos[0][2] = 0;
	
	pos[1][0] = 0;
	pos[1][1] = 1;
	pos[1][2] = 0;
	
	pos[2][0] = 1;
	pos[2][1] = 1;
	pos[2][2] = 0;
	
	pos[3][0] = 1;
	pos[3][1] = 2;
	pos[3][2] = 0;

	pos[4][0] = 1;
	pos[4][1] = 1;
	pos[4][2] = 1;

	pos[5][0] = 1;
	pos[5][1] = 0;
	pos[5][2] = 0;

	TorsionalCV* tortest = new TorsionalCV(1, 2, 3, 4, true);
	TorsionalCV* tortest2 = new TorsionalCV(1, 2, 3, 5, true);
	TorsionalCV* tortest3 = new TorsionalCV(1, 2, 3, 6, true);

	tortest->Initialize(snap);
	tortest2->Initialize(snap);
	tortest3->Initialize(snap);

	tortest->Evaluate(snap);
	tortest2->Evaluate(snap);
	tortest3->Evaluate(snap);


	auto& grad1 = tortest->GetGradient();
	EXPECT_TRUE(grad1[0].isApprox(Vector3{0, 0, 1}));
	EXPECT_TRUE(grad1[1].isApprox(Vector3{0, 0, -1}));
	EXPECT_TRUE(grad1[2].isApprox(Vector3{0, 0, -1}));
	EXPECT_TRUE(grad1[3].isApprox(Vector3{0, 0, 1}));
	EXPECT_TRUE(grad1[4].isApprox(Vector3{0, 0, 0}));
	EXPECT_TRUE(grad1[5].isApprox(Vector3{0, 0, 0}));

	auto& grad2 = tortest2->GetGradient();
	EXPECT_TRUE(grad2[0].isApprox(Vector3{0, 0, 1}));
	EXPECT_TRUE(grad2[1].isApprox(Vector3{0, 0, -1}));
	EXPECT_TRUE(grad2[2].isApprox(Vector3{0, 1, 0}));
	EXPECT_TRUE(grad2[3].isApprox(Vector3{0, 0, 0}));
	EXPECT_TRUE(grad2[4].isApprox(Vector3{0, -1, 0}));
	EXPECT_TRUE(grad2[5].isApprox(Vector3{0, 0, 0}));
	
	auto& grad3 = tortest3->GetGradient();
	EXPECT_TRUE(grad3[0].isApprox(Vector3{0, 0, 1}));
	EXPECT_TRUE(grad3[1].isApprox(Vector3{0, 0, -1}));
	EXPECT_TRUE(grad3[2].isApprox(Vector3{0, 0, 1}));
	EXPECT_TRUE(grad3[3].isApprox(Vector3{0, 0, 0}));
	EXPECT_TRUE(grad3[4].isApprox(Vector3{0, 0, 0}));
	EXPECT_TRUE(grad3[5].isApprox(Vector3{0, 0, -1}));

	EXPECT_NEAR(tortest->GetValue(), 3.14159, 0.01);
	EXPECT_NEAR(tortest2->GetValue(), -1.570796, 0.01);
	EXPECT_NEAR(tortest3->GetValue(), 0, 0.01);
}
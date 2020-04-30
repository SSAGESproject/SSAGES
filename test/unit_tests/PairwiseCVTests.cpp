#include "Tests.h"

#include "CVs/PairwiseCV.h"

using namespace SSAGES;

class PairwiseCVTests : public ::testing::Test{
protected:
	virtual void SetUp()
	{
		// cv1 = new PairwiseCV({1, 2, 3, 4},   new GaussianPK(1.0, 1.0));
		cv2 = new PairwiseCV({1, 2, 3}, {4}, new GaussianPK(2.0, 1.0));
		cv3 = new PairwiseCV({4}, {1, 2, 3}, new GaussianPK(2.0, 1.0));
		cv4 = new PairwiseCV({1, 2, 3}, {4}, new RationalSwitchPK(0, 1.5, 12, 24));
		cv5 = new PairwiseCV({4}, {1, 2, 3}, new RationalSwitchPK(0, 1.5, 12, 24));

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
		// delete cv1;
		delete cv2;
		delete cv3;
	}

	Snapshot* snapshot;
	// NearestNeighborsCV *cv1, *cv2, *cv3, *cv4, *cv5;
	PairwiseCV *cv2, *cv3, *cv4, *cv5;
	mxx::comm comm;
};

TEST_F(PairwiseCVTests, DefaultBehavior)
{
	// Gradient should sum to zero.
	double sum = 0;
/* Commented out until double counting issue resolved
 *  // Gaussian function tests
 *
 *	// Check CV1 calculations.
 *	cv1->Initialize(*snapshot);
 *	cv1->Evaluate(*snapshot);
 *	
 *	EXPECT_NEAR(cv1->GetValue(), 5.300240013212994, eps);
 *
 *	auto& grad1 = cv1->GetGradient(); 
 *	// NOTE: Gradient checks are only valid for single processor
 *	// but value checks for 2. This needs to be adjusted later.
 *	if(comm.size() == 1)
 *	{
 *		EXPECT_NEAR(grad1[0][0],  0.751412029625549, eps);
 *		EXPECT_NEAR(grad1[0][1], -0.641890504603067, eps);
 *		EXPECT_NEAR(grad1[0][2],  1.069998040552683, eps);
 *
 *		EXPECT_NEAR(grad1[1][0], -0.268240611427197, eps);
 *		EXPECT_NEAR(grad1[1][1],  0.497513059206846, eps);
 *		EXPECT_NEAR(grad1[1][2], -0.368931410427070, eps);
 *
 *		EXPECT_NEAR(grad1[2][0], -0.358278054625905, eps);
 *		EXPECT_NEAR(grad1[2][1], -0.178927048279713, eps);
 *		EXPECT_NEAR(grad1[2][2], -0.179351006346193, eps);
 *
 *		EXPECT_NEAR(grad1[3][0], -0.124893363572446, eps);
 *		EXPECT_NEAR(grad1[3][1],  0.323304493675933, eps);
 *		EXPECT_NEAR(grad1[3][2], -0.521715623779420, eps);
 *	}
 *
 *	for(auto& v : grad1)
 *		sum += v.array().sum();
 *	sum = mxx::allreduce(sum, std::plus<double>(), comm);
 *	EXPECT_NEAR(sum, 0, eps);
 */

	// Check CV2 calculations.
	cv2->Initialize(*snapshot);
	cv2->Evaluate(*snapshot);
	
	EXPECT_NEAR(cv2->GetValue(), 2.004802380585009, eps);

	auto& grad2 = cv2->GetGradient();
	if(comm.size() == 1)
	{
		EXPECT_NEAR(grad2[0][0], -0.149245537579517, eps);
		EXPECT_NEAR(grad2[0][1],  0.149245537579517, eps);
		EXPECT_NEAR(grad2[0][2], -0.149245537579517, eps);

		EXPECT_NEAR(grad2[1][0], -0.396342114343994, eps);
		EXPECT_NEAR(grad2[1][1],  0.000000000000000, eps);
		EXPECT_NEAR(grad2[1][2],  0.396342114343994, eps);

		EXPECT_NEAR(grad2[2][0],  0.000000000000000, eps);
		EXPECT_NEAR(grad2[2][1],  0.606530659712633, eps);
		EXPECT_NEAR(grad2[2][2],  0.000000000000000, eps);

		EXPECT_NEAR(grad2[3][0],  0.545587651923512, eps);
		EXPECT_NEAR(grad2[3][1], -0.755776197292151, eps);
		EXPECT_NEAR(grad2[3][2], -0.247096576764477, eps);
	}
	
	sum = 0.; 
	for(auto& v : grad2)
		sum += v.array().sum();
	sum = mxx::allreduce(sum, std::plus<double>(), comm);
	EXPECT_NEAR(sum, 0, eps);
	
	// Check CV3 calculations. CV3 should match CV2.
	cv3->Initialize(*snapshot);
	cv3->Evaluate(*snapshot);

	EXPECT_NEAR(cv3->GetValue(), cv2->GetValue(), eps);

	auto& grad3 = cv3->GetGradient();
	sum = 0.; 
	for(auto& v : grad3)
		sum += v.array().sum();
	sum = mxx::allreduce(sum, std::plus<double>(), comm);
	EXPECT_NEAR(sum, 0, eps);
	
	// Rational switching function tests
	
	cv4->Initialize(*snapshot);
    cv4->Evaluate(*snapshot);
    EXPECT_NEAR(cv4->GetValue(), 2.143319272335021, eps);

    auto& grad4 = cv4->GetGradient(); 
    // NOTE: Gradient checks are only valid for single processor
    // but value checks for 2. This needs to be adjusted later. 
    if(comm.size() == 1)
    {
        EXPECT_NEAR(grad4[0][0], 0.5130418964267503, eps);
        EXPECT_NEAR(grad4[0][1], -0.5130418964267503, eps);
        EXPECT_NEAR(grad4[0][2], 0.5130418964267503, eps);

        EXPECT_NEAR(grad4[1][0], 0.0014447794902723533, eps);
        EXPECT_NEAR(grad4[1][1], 0, eps);
        EXPECT_NEAR(grad4[1][2], -0.0014447794902723533, eps);

        EXPECT_NEAR(grad4[2][0], 0, eps);
        EXPECT_NEAR(grad4[2][1], -0.09107879745469, eps);
        EXPECT_NEAR(grad4[2][2], 0, eps);

        EXPECT_NEAR(grad4[3][0], -0.51448667591702268, eps);
        EXPECT_NEAR(grad4[3][1], 0.60412069388144074, eps);
        EXPECT_NEAR(grad4[3][2], -0.51159711693647802, eps);
    }

    // Gradient should sum to zero.
    sum = 0.; 
    for(auto& v : grad4)
        sum += v.array().sum();
    sum = mxx::allreduce(sum, std::plus<double>(), comm);
    EXPECT_NEAR(sum, 0, eps);

    // CV2 should match CV1.
    cv5->Initialize(*snapshot);
    cv5->Evaluate(*snapshot);

    EXPECT_NEAR(cv5->GetValue(), cv4->GetValue(), eps);

    auto& grad5 = cv5->GetGradient();
    sum = 0.; 
    for(auto& v : grad5)
        sum += v.array().sum();
    sum = mxx::allreduce(sum, std::plus<double>(), comm);
    EXPECT_NEAR(sum, 0, eps);

    if(comm.size() == 1)
    {
        EXPECT_NEAR(grad5[0][0], grad4[0][0], eps);
        EXPECT_NEAR(grad5[0][1], grad4[0][1], eps);
        EXPECT_NEAR(grad5[0][2], grad4[0][2], eps);

        EXPECT_NEAR(grad5[1][0], grad4[1][0], eps);
        EXPECT_NEAR(grad5[1][1], grad4[1][1], eps);
        EXPECT_NEAR(grad5[1][2], grad4[1][2], eps);

        EXPECT_NEAR(grad5[2][0], grad4[2][0], eps);
        EXPECT_NEAR(grad5[2][1], grad4[2][1], eps);
        EXPECT_NEAR(grad5[2][2], grad4[2][2], eps);

        EXPECT_NEAR(grad5[3][0], grad4[3][0], eps);
        EXPECT_NEAR(grad5[3][1], grad4[3][1], eps);
        EXPECT_NEAR(grad5[3][2], grad4[3][2], eps);
    }
}

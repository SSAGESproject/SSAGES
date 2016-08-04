#include "../src/Grids/Grid.h"
#include "../src/Grids/Grid1D.h"
#include "../src/Grids/Grid2D.h"
#include "../src/Grids/Grid3D.h"
#include "../src/Utility/BuildException.h"
#include "gtest/gtest.h"
#include <boost/mpi.hpp>

using namespace SSAGES;

// Test calculation up to accuracy of 10^-10
const double eps = 0.0000000001;

class GridsTest : public ::testing::Test {

public:
	virtual void SetUp() {

		Grid1n = new Grid1D({-1.2}, {0.8}, {false}, {11});
		Grid1p = new Grid1D({-1.2}, {0.8}, {true}, {11});

		Grid2nn = new Grid2D({-109.5, 50.1}, {-60, 100.1}, {false, false}, {11, 101});
		Grid2np = new Grid2D({-109.5, 50.1}, {-60, 100.1}, {false, true}, {11, 101});
		Grid2pn = new Grid2D({-109.5, 50.1}, {-60, 100.1}, {true, false}, {11, 101});
		Grid2pp = new Grid2D({-109.5, 50.1}, {-60, 100.1}, {true, true}, {11, 101});
	}
	
	virtual void TearDown() {
		delete Grid1n;
		delete Grid1p;

		delete Grid2pn;
		delete Grid2nn;
		delete Grid2np;
		delete Grid2pp;

	}

	Grid1D *Grid1n;
	Grid1D *Grid1p;

	Grid2D *Grid2pn;
	Grid2D *Grid2nn;
	Grid2D *Grid2np;
	Grid2D *Grid2pp;
};

TEST_F(GridsTest, GetSetGo) {

	double value;
	//1D Grid Grid1n = new Grid1D({-1.2}, {0.8}, {false}, {11});
	Grid1n->SetValue({2}, 5.2);
	value = Grid1n->GetValue({2});
	EXPECT_NEAR(value, 5.2, eps);
	EXPECT_NEAR(-1.0, Grid1n->GetLocation({1})[0], eps);

	//2D grid Grid2nn = new Grid2D({-109.5, 50.1}, {-60, 100.1}, {false, false}, {11, 101});
	Grid2nn->SetValue({2,3}, 4.5);
	value = Grid2nn->GetValue({2,3});
	EXPECT_NEAR(value, 4.5, eps);
	EXPECT_NEAR(-104.55, Grid2nn->GetLocation({1,1})[0], eps);
	EXPECT_NEAR(50.6, Grid2nn->GetLocation({1,1})[1], eps);

}

TEST_F(GridsTest, GetParameters) {

	std::vector<int> intsv;
	std::vector<double> doubsv;
	std::vector<bool> boolsv;

	int dim;

	//1D Grid
	doubsv = Grid1n->GetLower();
	EXPECT_DOUBLE_EQ(doubsv[0], -1.2);
	doubsv = Grid1n->GetUpper();
	EXPECT_DOUBLE_EQ(doubsv[0], 0.8);
	doubsv = Grid1n->GetSpacing();
	EXPECT_NEAR(doubsv[0], 0.20, eps);
	boolsv = Grid1n->GetPeriodic();
	EXPECT_FALSE(boolsv[0]);
	intsv = Grid1n->GetNumPoints();
	EXPECT_EQ(intsv[0], 11);
	dim = Grid1n->GetDimension();
	EXPECT_EQ(dim, 1);

	//2D Grid
	//Grid2nn = new Grid2D({-109.5, 50.1}, {-60, 100.1}, {false, false}, {11, 101});
	doubsv = Grid2nn->GetLower();
	EXPECT_DOUBLE_EQ(doubsv[0], -109.5);
	EXPECT_DOUBLE_EQ(doubsv[1], 50.1);
	doubsv = Grid2nn->GetUpper();
	EXPECT_DOUBLE_EQ(doubsv[0], -60);
	EXPECT_DOUBLE_EQ(doubsv[1], 100.1);
	boolsv = Grid2nn->GetPeriodic();
	EXPECT_FALSE(boolsv[0]);
	EXPECT_FALSE(boolsv[1]);
	intsv = Grid2nn->GetNumPoints();
	EXPECT_EQ(intsv[0], 11);
	EXPECT_EQ(intsv[1], 101);
	doubsv = Grid2nn->GetSpacing();
	EXPECT_NEAR(doubsv[0], 4.95, eps);
	EXPECT_NEAR(doubsv[1], 0.5, eps);
	dim = Grid2nn->GetDimension();
	EXPECT_EQ(dim, 2);
}

TEST_F(GridsTest, FlattenIndices) {
	EXPECT_EQ(0, FlattenIndices({0},Grid1n->GetNumPoints()));
	EXPECT_EQ(10, FlattenIndices({10},Grid1n->GetNumPoints()));
	EXPECT_EQ(3, FlattenIndices({3},Grid1n->GetNumPoints()));
	
	EXPECT_EQ(0, FlattenIndices({0,0},Grid2nn->GetNumPoints()));
	EXPECT_EQ(100, FlattenIndices({0,100},Grid2nn->GetNumPoints()));
	EXPECT_EQ(1010, FlattenIndices({10,0},Grid2nn->GetNumPoints()));
	EXPECT_EQ(1110, FlattenIndices({10,100},Grid2nn->GetNumPoints()));
	EXPECT_EQ(555, FlattenIndices({5,50},Grid2nn->GetNumPoints()));

}

TEST_F(GridsTest, OutOfBounds) {

	EXPECT_THROW(FlattenIndices({11},{10}),std::out_of_range);
	EXPECT_THROW(FlattenIndices({1,11},{10,10}),std::out_of_range);

	// Grid1n = new Grid1D({-1.2}, {0.8}, {false}, {11});
	EXPECT_EQ(10, Grid1n->GetIndices({1.0})[0]);
	EXPECT_EQ(0, Grid1n->GetIndices({-10.0})[0]);

	// Grid1p = new Grid1D({-1.2}, {0.8}, {true}, {11});
	EXPECT_EQ(1, Grid1p->GetIndices({1.15})[0]);
	EXPECT_EQ(7, Grid1p->GetIndices({-1.95})[0]);

	// Grid2nn = new Grid2D({-109.5, 50.1}, {-60, 100.1}, {false, false}, {11, 101});
	EXPECT_EQ(10, Grid2nn->GetIndices({100,100})[0]);
	EXPECT_EQ(100, Grid2nn->GetIndices({1000,1000})[1]);
	EXPECT_EQ(0, Grid2nn->GetIndices({-1000})[0]);
	EXPECT_EQ(0, Grid2nn->GetIndices({-10000})[1]);
	
	// Grid2pp = new Grid2D({-109.5, 50.1}, {-60, 100.1}, {true, true}, {11, 101});
	EXPECT_EQ(1, Grid2pp->GetIndices({-50.0,0})[0]); //4.95 spacing
	EXPECT_EQ(3, Grid2pp->GetIndices({0,102})[1]); //0.5 spacing
	EXPECT_EQ(9, Grid2pp->GetIndices({-117.5,0})[0]); //4.95 spacing
	EXPECT_EQ(99, Grid2pp->GetIndices({0,49.0})[1]); //0.5 spacing
}

int main(int argc, char *argv[])
{
    ::testing::InitGoogleTest(&argc, argv);
    int ret = RUN_ALL_TESTS();
    return ret;
}

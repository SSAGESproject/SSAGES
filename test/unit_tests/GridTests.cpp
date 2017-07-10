#include "gtest/gtest.h"
#include "Grids/Grid.h"
#include <random>

using namespace SSAGES;

class GridTest : public ::testing::Test {
protected:
    virtual void SetUp() 
    {
		grid = new Grid<double>({10, 10}, {-1.5, -2.0}, {10.0, 6.7}, {false, true});
	}

	virtual void TearDown() 
	{
		delete grid; 
    }

	Grid<double>* grid;
};

TEST_F(GridTest, ReadWriteGrid)
{
	// Fill grid with random numbers.
	std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(-1000., 1000.);

	auto* data1 = grid->data(); 
	auto size = grid->size();
	for(size_t i = 0; i < size; ++i)
		data1[i] = dis(gen);

	// Write grid to file.
	grid->WriteToFile("testgrid.dat");

	// Copy grid. 
	Grid<double> grid2 = *grid;
	
	// Clear old grid.
	for(size_t i = 0; i < size; ++i)
		data1[i] = 0;
	
	// Load grid back from file. 
	grid->LoadFromFile("testgrid.dat");

	// Compare to grid2.
	auto* data2 = grid2.data();

	for(size_t i = 0; i < size; ++i)
		EXPECT_NEAR(data1[i], data2[i], 1e-8);	
}

int main(int argc, char *argv[])
{
    ::testing::InitGoogleTest(&argc, argv);
    int ret = RUN_ALL_TESTS();
    return ret;
}

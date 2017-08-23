#include "gtest/gtest.h"
#include "Grids/Grid.h"
#include <random>

using namespace SSAGES;

class GridTest : public ::testing::Test {
protected:
    virtual void SetUp() 
    {
		grid = new Grid<double>({2, 2}, {-1.5, -2.0}, {10.0, 6.7}, {false, true});
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

TEST_F(GridTest, Interpolation)
{
	// Fill grid with random numbers.
	std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(-1000., 1000.);

	auto* data1 = grid->data(); 
	auto size = grid->size();
	for(size_t i = 0; i < size; ++i)
		data1[i] = dis(gen);

	// First, test a point within bounds. (3.0,3.0)
	double real_value = (grid->at({0,0})*std::abs(3-(-1.5+3*(10+1.5)/4))	*std::abs(3-(-2.0+3*(6.7+2.0)/4))	/ (((10+1.5)/2)*((6.7+2.0)/2))) + 
						(grid->at({1,0})*std::abs(3-(-1.5+(10+1.5)/4))	*std::abs(3-(-2.0+3*(6.7+2.0)/4))	/ (((10+1.5)/2)*((6.7+2.0)/2))) +
						(grid->at({0,1})*std::abs(3-(-1.5+3*(10+1.5)/4))	*std::abs(3-(-2.0+(6.7+2.0)/4))		/ (((10+1.5)/2)*((6.7+2.0)/2))) +
						(grid->at({1,1})*std::abs(3-(-1.5+(10+1.5)/4))	*std::abs(3-(-2.0+(6.7+2.0)/4))		/ (((10+1.5)/2)*((6.7+2.0)/2)));
						
	EXPECT_NEAR(real_value, grid->GetInterpolated({3.0,3.0}), 1e-8);
	
	// Next, test non-periodic index outside grid centers' span bounds. (8.0,3.0)
	real_value = 		(grid->at({1,0})*std::abs(3-(-2.0+3*(6.7+2.0)/4))	/ ((6.7+2.0)/2)) +
						(grid->at({1,1})*std::abs(3-(-2.0+(6.7+2.0)/4))		/ ((6.7+2.0)/2));
						
	EXPECT_NEAR(real_value, grid->GetInterpolated({8.0,3.0}), 1e-8);
	
	// Finally, test periodicity handling. (3.0,6.0)	
	real_value = 		(grid->at({0,0})*std::abs(3-(-1.5+3*(10+1.5)/4))	*std::abs(6-(-2.0+3*(6.7+2.0)/4))	/ (((10+1.5)/2)*((6.7+2.0)/2))) + 
						(grid->at({1,0})*std::abs(3-(-1.5+(10+1.5)/4))	*std::abs(6-(-2.0+3*(6.7+2.0)/4))	/ (((10+1.5)/2)*((6.7+2.0)/2))) +
						(grid->at({0,1})*std::abs(3-(-1.5+3*(10+1.5)/4))	*std::abs(6-(-2.0+5*(6.7+2.0)/4))		/ (((10+1.5)/2)*((6.7+2.0)/2))) +
						(grid->at({1,1})*std::abs(3-(-1.5+(10+1.5)/4))	*std::abs(6-(-2.0+5*(6.7+2.0)/4))		/ (((10+1.5)/2)*((6.7+2.0)/2)));
						
	EXPECT_NEAR(real_value, grid->GetInterpolated({3.0,6.0}), 1e-8);
}

int main(int argc, char *argv[])
{
    ::testing::InitGoogleTest(&argc, argv);
    int ret = RUN_ALL_TESTS();
    return ret;
}

#include "Tests.h"

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


TEST_F(GridTest, GridBounds)
{

	// Test GetLower, GetUpper with periodic and non-periodic bounds.
	EXPECT_NEAR(grid->GetLower(0),-1.5, 1e-8);	
	EXPECT_NEAR(grid->GetLower(1),-2.0, 1e-8);

	EXPECT_NEAR(grid->GetUpper(0),10.0, 1e-8);	
	EXPECT_NEAR(grid->GetUpper(1),6.7, 1e-8);

	// Test GetLower, GetUpper with periodic and non-periodic bounds with vector output.
	std::vector<double> lower = {-1.5,-2.0};
	std::vector<double> upper = {10.0,6.7};

	for(size_t i = 0; i < lower.size(); ++i)
	{
		EXPECT_NEAR(grid->GetLower()[i],lower[i], 1e-8);
		EXPECT_NEAR(grid->GetUpper()[i],upper[i], 1e-8);
	} 

}

TEST_F(GridTest, GridCenters)
{
	std::vector<double> bin00center = {1.375,-2.0};
	std::vector<double> bin01center = {1.375,6.7};
	std::vector<double> bin10center = {7.125,-2.0};
	std::vector<double> bin11center = {7.125,6.7};

	std::vector<double> output = grid->GetCoordinates({0,0});
	for(size_t i = 0; i < output.size(); ++i)
		EXPECT_NEAR(bin00center[i],output[i], 1e-8);

	output = grid->GetCoordinates({0,1});
	for(size_t i = 0; i < output.size(); ++i)
		EXPECT_NEAR(bin01center[i],output[i], 1e-8);

	output = grid->GetCoordinates({1,0});
	for(size_t i = 0; i < output.size(); ++i)
		EXPECT_NEAR(bin10center[i],output[i], 1e-8);

	output = grid->GetCoordinates({1,1});
	for(size_t i = 0; i < output.size(); ++i)
		EXPECT_NEAR(bin11center[i],output[i], 1e-8);
	
}

TEST_F(GridTest, SyncGrid)
{
	grid->at({-1.49,-1.9}) = 1.0;

	auto* data1 = grid->data();
	EXPECT_NEAR(data1[0],1.0, 1e-8);
	EXPECT_NEAR(data1[1],0.0, 1e-8);
	EXPECT_NEAR(data1[2],0.0, 1e-8);
	EXPECT_NEAR(data1[3],0.0, 1e-8);

	grid->syncGrid();

	EXPECT_NEAR(data1[0],1.0, 1e-8);
	EXPECT_NEAR(data1[1],0.0, 1e-8);
	EXPECT_NEAR(data1[2],1.0, 1e-8);
	EXPECT_NEAR(data1[3],0.0, 1e-8);

}



/*TEST_F(GridTest, Interpolation)
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
}*/

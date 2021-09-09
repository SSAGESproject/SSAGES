#include "Tests.h"

#include "CVs/CVManager.h"
#include "Methods/ABF.h"
#include "CVs/MockCV.h"

using namespace SSAGES;

class ABFTests : public ::testing::Test
{
protected:

    virtual void SetUp() 
	{

		// Set up dummy CVs.
		Vector3 grad1 = {1,1,0};

		CV1 = new MockCV(-0.5,{1,0,0},-2,2); // In bounds.
		CV2 = new MockCV(0.5,{0,1,0},-2,2); // In bounds.
		CV3 = new MockCV(-2,{0,0,1},-2,2); // Out of bounds.

		cvmanager.AddCV(CV1);
		cvmanager.AddCV(CV2);
		cvmanager.AddCV(CV3);

		// Set up three timepoints.
		snapshot1 = new Snapshot(comm, 0);
		snapshot2 = new Snapshot(comm, 0);
		snapshot3 = new Snapshot(comm, 0);

		// Unit conversion
		snapshot1->dMdt2f_ = 1;
		snapshot2->dMdt2f_ = 1;
		snapshot3->dMdt2f_ = 1;

		snapshot1->SetNumAtoms(1);
		snapshot2->SetNumAtoms(1);
		snapshot3->SetNumAtoms(1);

		auto& loc1 = snapshot1 ->GetPositions();
		auto& loc2 = snapshot2 ->GetPositions();
		auto& loc3 = snapshot3 ->GetPositions();

		loc1.resize(1);
		loc2.resize(1);
		loc3.resize(1);

		auto& vel1 = snapshot1 -> GetVelocities();
		auto& vel2 = snapshot2 -> GetVelocities();
		auto& vel3 = snapshot3 -> GetVelocities();

		vel1.resize(1);
		vel2.resize(1);
		vel3.resize(1);

		vel1[0][0] = 2.0;
		vel1[0][1] = 3.0;
		vel1[0][2] = 4.0;

		vel2[0][0] = 1.0;
		vel2[0][1] = 0.0;
		vel2[0][2] = -1.0;

		vel3[0][0] = 0.0;
		vel3[0][1] = -3.0;
		vel3[0][2] = -6.0;

		auto& forces1 = snapshot1 -> GetForces();
		auto& forces2 = snapshot2 -> GetForces();
		auto& forces3 = snapshot3 -> GetForces();

		forces1.resize(1);
		forces1[0].resize(3);
		forces1[0] = {0,0,0};

		forces2.resize(1);
		forces2[0].resize(3);
		forces2[0] = {0,0,0};

		forces3.resize(1);
		forces3[0].resize(3);
		forces3[0] = {0,0,0};

		auto& mass1 = snapshot1 -> GetMasses();
		auto& mass2 = snapshot2 -> GetMasses();
		auto& mass3 = snapshot3 -> GetMasses();

		mass1.resize(1);
		mass2.resize(1);
		mass3.resize(1);

		mass1[0] = 2.0;
		mass2[0] = 2.0;
		mass3[0] = 2.0;

		// unsigned int mpiid1 = snapshot1 -> GetWalkerID();
		// unsigned int mpiid2 = snapshot2 -> GetWalkerID();
		// unsigned int mpiid3 = snapshot3 -> GetWalkerID();

		// mpiid1 = 0;
		// mpiid2 = 0;
		// mpiid3 = 0;

		// Set up ABF method.
		std::vector<std::vector<double>> histdetails;
		std::vector<std::vector<double>> restraint;
		std::vector<bool> isperiodic;
		std::vector<std::vector<double>> periodicboundaries;
					
		for(size_t i=0; i < 3; ++i)
		{
			std::vector<double> temp1 = {-1, 1, 2};
			std::vector<double> temp2 = {-1.5, 1.5, 10};
			histdetails.push_back(temp1);
			restraint.push_back(temp2);
			isperiodic.push_back(false);
			periodicboundaries.push_back(std::vector<double>(2, 0.0));
		}
		
		Grid<int> *N;
		Grid<int> *Nworld;
		std::vector<Grid<double>*> F(3);
		std::vector<Grid<double>*> Fworld(3);
		
		N= new Grid<int>({2,2,2}, {-1,-1,-1}, {1,1,1}, isperiodic);
		Nworld= new Grid<int>({2,2,2}, {-1,-1,-1}, {1,1,1}, isperiodic);
			
		for(auto& grid : F)
		{
			grid= new Grid<double>({2,2,2}, {-1,-1,-1}, {1,1,1}, isperiodic);
		}
		for(auto& grid : Fworld)
		{
			grid= new Grid<double>({2,2,2}, {-1,-1,-1}, {1,1,1}, isperiodic);
		}

		Method = new ABF(world, // MPI global communicator
						comm,  // MPI local communicator
						N, 
						Nworld, 
						F, 
						Fworld, 
						restraint, // Min, max, spring constant
						isperiodic, // Periodicity in CVs
						periodicboundaries, // Position of periodic boundaries
						5, // min hists before applying bias
						false, // mass weighing
						"testout", // filename
						true, // overwrite testout
						"Nworld", //Nworld_filename
						"Fworld_cv",  //Fworld_filename
						histdetails, // hist details
						10, // Backup interval
						1, // timestep
						1); // frequency
					
	}

	virtual void TearDown() 
	{
		delete snapshot1;
		delete snapshot2;
		delete snapshot3;

		delete Method;
	}

	mxx::comm world;
	mxx::comm comm;

	Snapshot* snapshot1;
	Snapshot* snapshot2;
	Snapshot* snapshot3;

	MockCV* CV1;
	MockCV* CV2;
	MockCV* CV3;

	ABF* Method;	

	CVManager cvmanager;
};

TEST_F(ABFTests,Initialization)
{
	// Initialize CVs.	
	CV1->MockCV::Initialize(*snapshot1);
	CV2->MockCV::Initialize(*snapshot1);
	CV3->MockCV::Initialize(*snapshot1);

	// Initialize Method.
	Method->ABF::PreSimulation(snapshot1, cvmanager);

	// Test dimensionality is correct.
	int dim = 3;
	
	// Check for correct initialization.	
	EXPECT_EQ(Method->N_->size(), 8);
	EXPECT_EQ(Method->F_.size(), 3);

	/*
	for(size_t i = 0; i < dim ; ++i)
	{
		EXPECT_TRUE(Method->N_ == 0);

		EXPECT_TRUE((Method->F_)[0][i] == 0);
		EXPECT_TRUE((Method->F_)[1][i] == 0);
		EXPECT_TRUE((Method->F_)[2][i] == 0);
	}
	*/
	
}

/*
TEST_F(ABFTests,Getters_Setters)
{
	Method->ABF::PreSimulation(snapshot1, cvmanager);

	// Check Setters
	Eigen::VectorXd F(2);
	F << 1,1;
	F_[0] = F;
	EXPECT_TRUE(Method->N_[0] == 1);
	EXPECT_TRUE(Method->N_[1] == 1);
	EXPECT_TRUE(Method->F_[0] == 1);
	EXPECT_TRUE(Method->F_[1] == 1);
}
*/

TEST_F(ABFTests,CV_outofboundscheck)
{
	// Initialize CVs.	
	CV1->MockCV::Initialize(*snapshot1);
	CV2->MockCV::Initialize(*snapshot1);
	CV3->MockCV::Initialize(*snapshot1);

	// Initialize Method.
	Method->ABF::PreSimulation(snapshot1, cvmanager);
	int dim = 3;
	
	// Take one MD step.
	Method->ABF::PostIntegration(snapshot1, cvmanager);

	// Test whether bound checking works. There should be no hits as 3rd CV was out of bounds.
	
	for(auto& point : *(Method->N_))
	{
		EXPECT_EQ(point, 0);
	}
	/*
	{
		EXPECT_TRUE(Method->_N[i] == 0);

		EXPECT_TRUE(Method->_F[0][i] == 0);
		EXPECT_TRUE(Method->_F[1][i] == 0);
		EXPECT_TRUE(Method->_F[2][i] == 0);
	}
	*/

	// The only bias force should be from CV restraints.
	double bias = (CV3->GetValue()+1.5)*10;

	for(size_t i = 0; i < Method->biases_[0].size() ; ++i)
		EXPECT_NEAR(-Method->biases_[0][i], bias*CV3->GetGradient()[0][i], eps);
}

TEST_F(ABFTests,MD_steps_inboundscheck)
{
	// Initialize CVs.	
	CV1->MockCV::Initialize(*snapshot1);
	CV2->MockCV::Initialize(*snapshot1);
	CV3->MockCV::Initialize(*snapshot1);

	// Initialize Method.
	Method->ABF::PreSimulation(snapshot1, cvmanager);
	int dim = 3;	
	
	// Take an MD step (Out of bounds).
	Method->ABF::PostIntegration(snapshot1, cvmanager);

	// Put CVs in bound.
	CV3->val_ = -1;

	// Take two more MD steps (in bounds).
	Method->ABF::PostIntegration(snapshot2, cvmanager);
	Method->ABF::PostIntegration(snapshot3, cvmanager);

	// Check that hits are recorded correctly
	EXPECT_EQ(Method->N_->at({0,1,0}), 2);

	EXPECT_NEAR(Method->F_[0]->at({0,1,0}), -8.0, eps);
	EXPECT_NEAR(Method->F_[1]->at({0,1,0}), -20.4, eps);
	EXPECT_NEAR(Method->F_[2]->at({0,1,0}), -32.8, eps);

	// Check that the biases are accurate
	EXPECT_NEAR(Method->biases_[0][0], 1.6, eps);
	EXPECT_NEAR(Method->biases_[0][1], 4.08, eps);
	EXPECT_NEAR(Method->biases_[0][2], 6.56, eps);

	// Check that the biases added to forces correctly
	EXPECT_NEAR(snapshot3->GetForces()[0][0], 1.6, eps);
	EXPECT_NEAR(snapshot3->GetForces()[0][1], 4.08, eps);
	EXPECT_NEAR(snapshot3->GetForces()[0][2], 6.56, eps);
}

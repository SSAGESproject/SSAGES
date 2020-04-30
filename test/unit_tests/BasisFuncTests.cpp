#include "Tests.h"

#include "CVs/CVManager.h"
#include "Methods/BasisFunc.h"
#include "CVs/MockCV.h"

using namespace SSAGES;

class BasisFuncTests : public ::testing::Test {
protected:
    virtual void SetUp()
    {
		// Set up dummy CVs.
		Vector3 grad1 = {1,1,0};

		CV1 = new MockCV(-0.5,{1,0,0},-1.5,2); // In bounds.
		CV2 = new MockCV(0.5,{0,1,0},-2,1.5); // In bounds.

		cvmanager.AddCV(CV1);
		cvmanager.AddCV(CV2);

		// Set up three timepoints.
		snapshot1 = new Snapshot(comm, 0);
		snapshot2 = new Snapshot(comm, 0);
		snapshot3 = new Snapshot(comm, 0);		

		snapshot1->SetNumAtoms(1);
		snapshot2->SetNumAtoms(1);
		snapshot3->SetNumAtoms(1);

        // Set up the boltzmann constant for all the snapshots
		snapshot1->SetKb(1.0);
		snapshot2->SetKb(1.0);
		snapshot3->SetKb(1.0);

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

        h = new Grid<unsigned int>({10,10}, {-1.5,-2.0}, {2.0,1.5}, {true, false});
        b = new Grid<double>({10,10}, {-1.5,-2.0}, {2.0,1.5}, {true, false});
        f = new Grid<std::vector<double>> ({10,10}, {-1.5,-2.0}, {2.0,1.5}, {true, false});
        std::vector<BasisFunction*> functions;
        functions.push_back(new Legendre(10,10));
        functions.push_back(new Chebyshev(10,-2.0,1.5,10));
        std::vector<double> restr = {10,10};
        std::vector<double> bndUp = {2.1,1.6};
        std::vector<double> bndLw = {-1.6,-2.1};
        
        Method = new BFS(world, // MPI global communicator
                         comm, // MPI local communicator
                         h, f, b, // histogram of visited states, bias, force grids
                         functions, // basis function vector
                         restr, bndUp, bndLw, // boundary positions and restraints
                         3, // frequency of bias updates
                         1, // frequency
                         "testout", // filename
                         1.0, // input temperature for poorly defined systems
                         eps, // tolerance
                         1.0, // weighting factor
                         false); // convergence criteria
    }

    virtual void TearDown() 
    {
		delete snapshot1;
		delete snapshot2;
		delete snapshot3;

        delete h;
        delete f;
        delete b;

		delete Method;
    }

	mxx::comm world;
	mxx::comm comm;

	Snapshot* snapshot1;
	Snapshot* snapshot2;
	Snapshot* snapshot3;

    Grid<unsigned int>* h;
    Grid<double>* b;
    Grid<std::vector<double>>* f;

	MockCV* CV1;
	MockCV* CV2;

    BFS* Method;

	CVManager cvmanager;
};

TEST_F(BasisFuncTests,Initialization)
{
	// Initialize CVs.	
	CV1->MockCV::Initialize(*snapshot1);
	CV2->MockCV::Initialize(*snapshot1);

	// Initialize Method.
	Method->BFS::PreSimulation(snapshot1, cvmanager);

	// Test dimensionality is correct.
	int dim = 3;
	// Check for correct initialization.
	EXPECT_TRUE(Method->unbias_.size() == 110);
	EXPECT_TRUE(Method->f_->size() == 110);

    // Check to make sure all histograms have a zero initialization
    size_t i = 0;
	for(Grid<unsigned int>::iterator it = Method->h_->begin(); it != Method->h_->end(); ++it, ++i)
	{
		EXPECT_TRUE(Method->b_->at(it.coordinates()) == 0);
		EXPECT_TRUE(Method->h_->at(it.coordinates()) == 0);
	}

    // Check the evaluator was constructed correctly
    EXPECT_TRUE(Method->evaluator_.functions_.size() == 2);
    EXPECT_TRUE(Method->evaluator_.coeff_.size() == 121); //10 + 1 coeffs
    EXPECT_NEAR(Method->evaluator_.lookup_[0].values[73],0.224073,0.00001);
    EXPECT_NEAR(Method->evaluator_.lookup_[1].derivs[41],0.65380758017492635,0.000001);
}

TEST_F(BasisFuncTests,InBounds)
{
	// Initialize CVs.	
	CV1->MockCV::Initialize(*snapshot1);
	CV2->MockCV::Initialize(*snapshot1);

	// Initialize Method.
	Method->BFS::PreSimulation(snapshot1, cvmanager);

    // Take three steps
    Method->BFS::PostIntegration(snapshot1, cvmanager);

    size_t i = 0;

	// Check the coefficients as well
    EXPECT_NEAR(Method->unbias_[80],2.0,eps);
	EXPECT_NEAR(Method->evaluator_.coeff_[11].value,0.009171649,eps);
	EXPECT_NEAR(Method->evaluator_.coeff_[36].value,-0.0476908749576,eps);
}

TEST_F(BasisFuncTests,OutBounds)
{
	// Initialize CVs.	
	CV1->MockCV::Initialize(*snapshot1);
	CV2->MockCV::Initialize(*snapshot1);

	// Initialize Method.
	Method->BFS::PreSimulation(snapshot1, cvmanager);
    CV2->val_ = -3;
    CV1->val_ = -3.5;

    // Take a single step
    Method->BFS::PostIntegration(snapshot1, cvmanager);

    EXPECT_FALSE(Method->bounds_);

	for(Grid<unsigned int>::iterator it = Method->h_->begin(); it != Method->h_->end(); ++it)
	{
		EXPECT_TRUE(Method->h_->at(it.coordinates()) == 0);
	}
}

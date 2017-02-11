#include "../hooks/lammps/fix_ssages.h"
#include "gtest/gtest.h"
#include <boost/mpi.hpp>

using namespace SSAGES;

TEST(UCVector, ConvertToUCVector)
{
	boost::mpi::communicator comm;

	std::array<double,6> box;
	box[0] = 35.3;
	box[1] = 29.4;
	box[2] = 11.8;

	box[3] = box[4] = box[5] = 0.0;

	const auto& newbox = LAMMPS_NS::ConvertToLatticeConstant(box);

	EXPECT_NEAR(newbox[0], 35.3, 1E-5 );
	EXPECT_NEAR(newbox[1], 29.4, 1E-5);
	EXPECT_NEAR(newbox[2], 11.8, 1E-5);
	EXPECT_NEAR(newbox[3], 1.5708, 0.001);
	EXPECT_NEAR(newbox[4], 1.5708, 0.001);
	EXPECT_NEAR(newbox[5], 1.5708, 0.001);

	box[0] = 35.3;
	box[1] = 29.4;
	box[2] = 11.8;

	box[3] = 2;
	box[4] = 2;
	box[5] = 2;

	const auto& newbox2 = LAMMPS_NS::ConvertToLatticeConstant(box);

	//Need to finish function here
}

// TEST(UCVector, GatherLAMMPSVectors)
// {
	// boost::mpi::communicator comm;

	// Hook* hook_;

	// std::shared_ptr<LAMMPS> lammps_;

	// // Silence of the lammps.
	// char **largs = (char**) malloc(sizeof(char*) * 5);
	// for(int i = 0; i < 5; ++i)
	// 	largs[i] = (char*) malloc(sizeof(char) * 1024);
	// sprintf(largs[0], " ");
	// sprintf(largs[1], "-screen");
	// sprintf(largs[2], "none");
	// sprintf(largs[3], "-log");
	// sprintf(largs[4], "none");

	// lammps_ = std::make_shared<LAMMPS>(5, largs, MPI_Comm(comm));

	// lammps_->input->one(("region simple block 0.5 35.8 0.6 30.0 -1.2 10.6").c_str());

	// if(!(hook_ = dynamic_cast<Hook*>(lammps_->modify->fix[fid])))
	// {
	// 	throw BuildException({"Unable to dynamic cast hook on node " + std::to_string(world_.rank())});			
	// }

	//Need to finish Unit test here
// }
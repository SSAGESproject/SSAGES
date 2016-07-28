#pragma once

#include "UnitCellConversion.h"


namespace SSAGES
{
//! Helper function to give the shortest distance between two points, accounting for PBCs.
	/*!
	 *  Uses the fractional coordinates - converts to fractional, then for each dimension accounts for PBCs (Xfrac < 0.5)
	 *  then converts back to cartesian. Returns a Vector3.
	 *
	 *  \ingroup Utility
	 */

	Vector3 NearestNeighbor(const std::array<double, 6>& LatticeConstants, const Vector3& a, const Vector3& b)
	{
		Vector3 parallelogram{{LatticeConstants[0],LatticeConstants[1],LatticeConstants[2]}};
		Vector3 angles{{LatticeConstants[3],LatticeConstants[4],LatticeConstants[5]}};

		std::vector<std::vector<double>> c2f(3,std::vector<double>(3));
		std::vector<std::vector<double>> f2c(3,std::vector<double>(3));

		c2f = uc_c2f(parallelogram, angles);
		f2c = uc_f2c(parallelogram, angles);

		// Convert the cartesian distance to fractional coordinates
		Vector3 del{{	 c2f[0][0]*(a[0]-b[0]) + c2f[0][1]*(a[1]-b[1]) + c2f[0][2]*(a[2]-b[2])
				,c2f[1][1]*(a[1]-b[1]) + c2f[1][2]*(a[2]-b[2])
				,c2f[2][2]*(a[2]-b[2]) }};
	
		// Take shortest distance through PBCs
			for(size_t i = 0; i<del.size(); ++i)
				{
				while(del[i] < -0.5){
					del[i] += 1.0;
					}
				while(del[i] > 0.5){
					del[i] -= 1.0;
					}
				}

		// Convert back to cartesian coordinates
			del = {{ f2c[0][0]*del[0] + f2c[0][1]*del[1] + f2c[0][2]*del[2]
				,f2c[1][1]*del[1] + f2c[1][2]*del[2]
				,f2c[2][2]*del[2]}};

		return del;
	}

}

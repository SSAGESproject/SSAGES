#pragma once

#include "UnitCellConversion.h"


namespace SSAGES
{
//! Helper function to unwrap coordinates .
	/*!
	 *  Converts Cartesian to fractional, adjusts the fractional coordinates to appropriate images,
	 *  then converts back to Cartesian. Returns unwrapped coordinates; Vector 3. 
	 *  \ingroup Utility
	 */

	Vector3 UnwrapCoordinates(const std::array<double, 6>& LatticeConstants, const Vector3& image, const Vector3& coord)

		//double unwrap_coord(const Vector3& coord, Vector3& image, Vector3& b, Vector3& a)
		{
			
			Vector3 parallelogram{{LatticeConstants[0],LatticeConstants[1],LatticeConstants[2]}};
			Vector3 angles{{LatticeConstants[3],LatticeConstants[4],LatticeConstants[5]}};

			std::vector<std::vector<double>> c2f(3,std::vector<double>(3));
			std::vector<std::vector<double>> f2c(3,std::vector<double>(3));

			c2f = uc_c2f(parallelogram, angles);
			f2c = uc_f2c(parallelogram, angles);

		// Convert Cartesian to fractional coordinates
			Vector3 fract_coord{{c2f[0][0]*(coord[0]]) + c2f[0][1]*(coord[1]) + c2f[0][2]*(coord[2])
							,c2f[1][1]*(coord[1]) + c2f[1][2]*(coord[2])
							,c2f[2][2]*(coord[2]) }};

		// Add pertinent images

			Vector3 un_fract_coord{{fract_coord[0] + image[0],
									fract_coord[1] + image[1],
									fract_coord[2] + image[2]}};
		// Convert back to cartesian coordinates
			Vector3 un_coord{{f2c[0][0]*un_fract_coord[0] + f2c[0][1]*un_fract_coord[1] + f2c[0][2]*un_fract_coord[2]
						,f2c[1][1]*un_fract_coord[1] + f2c[1][2]*un_fract_coord[2]
						,f2c[2][2]*un_fract_coord[2]}};

			return un_coord


		}
}

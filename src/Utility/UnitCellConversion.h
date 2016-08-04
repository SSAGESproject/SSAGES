#pragma once

#include <vector>
#include <array>
#include <math.h>

namespace SSAGES
{
	//! Helper function to find matrix to convert to fractional coordinates
	//! from cartesian.
	/*!
	 * \param LC Lattice constants of simulation box.
	 *
	 * \return Matrix to convert cartesian to fractional coordinates.
	 *
	 * c2f converts from cartesian to fractional coordinates, given lattice
	 * parameters a, b, c and alpha, beta and gamma.
	 *
	 * \ingroup Utility
	 */
		std::array<std::array<double,3>,3> uc_c2f(const std::array<double, 6>& LC)
		{
			std::array<std::array<double,3>,3> c2f;
			double v=0;

			v = sqrt(1-cos(LC[3])*cos(LC[3])-cos(LC[4])*cos(LC[4])-cos(LC[5])*cos(LC[5])+2*cos(LC[3])*cos(LC[4])*cos(LC[5]));

			c2f[0][0] = 1/LC[0];
			c2f[1][0] = 0;
			c2f[2][0] = 0;
			c2f[0][1] = -c2f[0][0]*cos(LC[5])/sin(LC[5]);
			c2f[1][1] = 1/(LC[1]*sin(LC[5]));
			c2f[2][1] = 0;
			c2f[0][2] = (cos(LC[3])*cos(LC[5]) - cos(LC[4]))/(LC[0]*v*sin(LC[5]));
			c2f[1][2] = (cos(LC[4])*cos(LC[5]) - cos(LC[3]))/(LC[1]*v*sin(LC[5]));
			c2f[2][2] = sin(LC[5])/(LC[2]*v);

			return c2f;

		}
	
		//! Helper function to find matrix to convert to cartesian coordinates
		//! from fractional.
		/*!
		 * \param LC Lattice constants (a, b, c, alpha, beta, gamma).
		 *
		 * \return Matrix for conversions from fractional to cartesian coordinates.
		 *
		 * \ingroup Utility
		 */
		std::array<std::array<double,3>,3> uc_f2c(const std::array<double, 6>& LC)
		{
			std::array<std::array<double,3>,3> f2c;

			f2c[0][0] = LC[0];
			f2c[1][0] = 0;
			f2c[2][0] = 0;
			f2c[0][1] = LC[1]*cos(LC[5]);
			f2c[1][1] = LC[1]*sin(LC[5]);
			f2c[2][1] = 0;
			f2c[0][2] = LC[2]*cos(LC[4]);
			f2c[1][2] = (LC[2]*LC[1]*cos(LC[3]) - f2c[0][1] * f2c[0][2])/ f2c[1][1];
			f2c[2][2] = sqrt(LC[2]*LC[2] - f2c[0][2]*f2c[0][2] - f2c[1][2]*f2c[1][2]);

			return f2c;

		}
	
}

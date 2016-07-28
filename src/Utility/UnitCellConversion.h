#pragma once

#include <vector>
#include <array>
#include <math.h>

namespace SSAGES
{

//! Helper function to convert to and from fractional coordinates using lattice parameters.
	/*!
	 *  c2f converts from cartesian to fractional coordinates, given lattice parameters a, b, c and alpha, beta and gamma.
	 *  f2c converts from fractional to cartesian coordinates, given lattice parameters a, b, c and alpha, beta and gamma.
	 *  Each function returns a Vector3.
	 *  \ingroup Utility
	 */


//! Helper function to find matrix to convert to fractional coordinates from cartesian.
		std::vector<std::vector<double> > uc_c2f(const Vector3& b, const Vector3& a)
		{
			std::vector<std::vector<double>> c2f(3,std::vector<double>(3));
			double ax=0;
			double ay=0;
			double az=0;
			double bx=0;
			double by=0;
			double bz=0; 
			double cx=0;
			double cy=0;
			double cz=0;
			double v=0;

			v = sqrt(1-cos(a[0])*cos(a[0])-cos(a[1])*cos(a[1])-cos(a[2])*cos(a[2])+2*cos(a[0])*cos(a[1])*cos(a[2]));

			ax = 1/b[0];
			ay = 0;
			az = 0;
			bx = -ax*cos(a[2])/sin(a[2]);
			by = 1/(b[1]*sin(a[2]));
			bz = 0;
			cx = (cos(a[0])*cos(a[2]) - cos(a[1]))/(b[0]*v*sin(a[2]));
			cy = (cos(a[1])*cos(a[2]) - cos(a[0]))/(b[1]*v*sin(a[2]));
			cz = sin(a[2])/(b[2]*v);

			c2f[0][0] = ax;
			c2f[1][0] = ay;
			c2f[2][0] = az;
			c2f[0][1] = bx;
			c2f[1][1] = by;
			c2f[2][1] = bz;
			c2f[0][2] = cx;
			c2f[1][2] = cy;
			c2f[2][2] = cz;

			return c2f;

		}
	
		//! Helper function to find matrix to convert to cartesian coordinates from fractional.
		std::vector< std::vector<double> > uc_f2c(const Vector3& b, const Vector3& a)
		{
			std::vector<std::vector<double>> f2c(3,std::vector<double>(3));
			double ax;
			double ay;
			double az;
			double bx;
			double by;
			double bz; 
			double cx;
			double cy;
			double cz;

			ax = b[0];
			ay = 0;
			az = 0;
			bx = b[1]*cos(a[2]);
			by = b[1]*sin(a[2]);
			bz = 0;
			cx = b[2]*cos(a[1]);
			cy = (b[2]*b[1]*cos(a[0]) - bx * cx)/ by;
			cz = sqrt(b[2]*b[2] - cx*cx - cy*cy);

			f2c[0][0] = ax;
			f2c[1][0] = ay;
			f2c[2][0] = az;
			f2c[0][1] = bx;
			f2c[1][1] = by;
			f2c[2][1] = bz;
			f2c[0][2] = cx;
			f2c[1][2] = cy;
			f2c[2][2] = cz;

			return f2c;

		}
	
}

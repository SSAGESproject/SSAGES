#pragma once

#include <vector>
#include <array>
#include <math.h>

namespace SSAGES
{

//! Helper functions to carry out vector operations.
	/*!
	 *  Has functions to return norm, norm squared, dot product and cross product of vectors.
	 *  \ingroup Utility
	 */
		
		//! Helper function to compute the norm of a vector.
		/*!
		 * \param v Three dimensional vector.
		 * \return Norm of the vector.
		 */
		double norm(const Vector3& v)
		{
			return sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
		}

		//! Helper function to compute the norm squared of a vector.
		/*!
		 * \param v Three dimensional vector.
		 * \return Norm squared of the vector.
		 */
		double norm2(const Vector3& v)
		{
			return (v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
		}

		//! Helper function to compute the dot product of two vectors.
		/*!
		 * \param v Three dimensional vector.
		 * \param w Three dimensional vector.
		 * \return Dot product v.w.
		 */
		double DotProduct(const Vector3& v, const Vector3& w)
		{
			return (v[0]*w[0] + v[1]*w[1] + v[2]*w[2]);
		}

		//! Helper function to compute the cross product of two vectors.
		/*!
		 * \param v Three dimensional vector.
		 * \param w Three dimensional vector.
		 * \return Dot product v x w.
		 */
		Vector3 CrossProduct(const Vector3& v, const Vector3& w)
		{
			Vector3 Cross;
			Cross[0] = v[1]*w[2] - v[2]*w[1];
			Cross[1] = v[2]*w[0] - v[0]*w[2];
			Cross[2] = v[0]*w[1] - v[1]*w[0];

			return Cross;
		}

		//! Helper function to compute the outer product of two vectors.
		/*!
		 * \param v Three dimensional vector.
		 * \param w Three dimensional vector.
		 * \return Outer product v x w.
		 */
		std::array<std::array<double, 3>,3> TensorProduct(const Vector3& v, const Vector3& w)
		{
			std::array<std::array<double, 3>,3> Tensor;
			for(size_t i; i < Tensor.size(); ++i)
				for(size_t j; j < Tensor[i].size(); ++j)
					Tensor[i][j] = v[i]*w[j];
			return Tensor;
		}
		
		//! Helper function to compute the determinant |ABC| of three vectors.
		/*!
		 * \param v Three dimensional vector.
		 * \param w Three dimensional vector.
		 * \param z Three dimensional vector.
		 * \return Determinant of matrix [v,w,z] = |vwz|.
		 */
		double Determinant(const Vector3& v, const Vector3& w , const Vector3& z)
		{
			return (v[0]*(w[1]*z[2]-z[1]*w[2]) -w[0]*(v[1]*z[2]-z[1]*v[2]) +z[0]*(v[1]*w[2]-w[1]*v[2]));
		}
		



	
}

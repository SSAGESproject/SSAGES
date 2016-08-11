#pragma once

#include "UnitCellConversion.h"
#include "types.h"

namespace SSAGES
{
    //! Helper function to unwrap coordinates.
    /*!
     * \param LatticeConstants Lattice Constants (a, b, c, alpha, beta, gamma).
     * \param coord Coordinates.
     * \param image Mirror images in the three dimensions.
     *
     * \return Vector3 containing the unwrapped coordinates.
     *
     * This function takes a set of (wrapped) coordinates and the information in
     * which mirror box the atom resides in. It returns the unwrapped atom
     * position. This is achieved by converting the Cartesian coordinates to
     * fractional coordinates, adjusting the fractional coordinates to the
     * appropriate mirror box, and then converting back to Cartesian coordinates.
     *
     * \note This function does not require the initial coordinates to be within
     *       the simulation box.
     *
     * \note This function assumes that the simulation box is periodic in all
     *       three dimensions.
     *
     * \note The values for the mirror images can take fractional values. The
     *       use of fractional images is highly discouraged and will be made
     *       invalid soon.
     *
     * \ingroup Utility
     */
    Vector3 UnwrapCoordinates(const std::array<double, 6>& LatticeConstants,
                              const Vector3& coord,
                              const Vector3& image)
    {
        auto c2f = uc_c2f(LatticeConstants);
        auto f2c = uc_f2c(LatticeConstants);

        // Convert Cartesian to fractional coordinates
        Vector3 fract(
            c2f[0][0]*(coord[0]) + c2f[0][1]*(coord[1]) + c2f[0][2]*(coord[2]),
            c2f[1][1]*(coord[1]) + c2f[1][2]*(coord[2]),
            c2f[2][2]*(coord[2])
        );

        // Add pertinent images
        Vector3 un_fract(fract[0] + image[0],
                         fract[1] + image[1],
                         fract[2] + image[2]);

        // Convert back to cartesian coordinates
        Vector3 un_coord(
            f2c[0][0]*un_fract[0] + f2c[0][1]*un_fract[1] + f2c[0][2]*un_fract[2],
            f2c[1][1]*un_fract[1] + f2c[1][2]*un_fract[2],
            f2c[2][2]*un_fract[2]
        );

        return un_coord;
    }
}

/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef DEFAULTPOTENTIAL_HPP_
#define DEFAULTPOTENTIAL_HPP_

#include "Tensor.hpp"
#include "simd.hpp"

class DefaultPotential
{
  public:
    //! The constructor
    DefaultPotential() {}

    //! Set the potential function for the scalar field here to zero
    template <class data_t, template <typename> class vars_t>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
    compute_potential(data_t &V_of_phi, data_t &dVdphi,
                      const vars_t<data_t> &vars) const
    {
        // The potential value at phi
        V_of_phi = 0.0;

        // The potential gradient at phi
        dVdphi = 0.0;
    }
};

#endif /* DEFAULTPOTENTIAL_HPP_ */

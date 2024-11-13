#include "StateVariables.hpp"
#include <AMReX_FArrayBox.H>
#include <AMReX_MultiFab.H>

#include "simd.hpp"
#include <cmath>

#ifndef POTENTIAL_HPP_
#define POTENTIAL_HPP_

class Potential
{
  private:
    amrex::Real m_mass;

  public:
    Potential(const amrex::Real mass) : m_mass(mass){};
    virtual ~Potential() = default;

    template <class data_t, template <typename> class vars_t>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE static void
    compute_phi_sq(data_t &V_of_phi, data_t &dVdphi,
                   const vars_t<data_t> &vars) 
    {

        V_of_phi = 0.5 * m_mass * m_mass * vars.phi * vars.phi;

        dVdphi = m_mass * m_mass * vars.phi;
    }

    template <class data_t, template <typename> class vars_t>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE static void
    compute_sine_gordon(data_t &V_of_phi, data_t &dVdphi,
                        const vars_t<data_t> &vars)
    {

        V_of_phi = std::sin(vars.phi);

        dVdphi = std::cos(vars.phi);
    }
};

#endif /* POTENTIAL_HPP_ */

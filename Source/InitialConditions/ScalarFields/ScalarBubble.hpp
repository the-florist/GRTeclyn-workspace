/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef SCALARBUBBLE_HPP_
#define SCALARBUBBLE_HPP_

#include "Cell.hpp"
#include "Coordinates.hpp"
#include "MatterCCZ4RHS.hpp"
#include "ScalarField.hpp"
#include "StateVariables.hpp" //This files needs NUM_VARS - total no. components
#include "Tensor.hpp"
#include "VarsTools.hpp"
#include "simd.hpp"

//! Class which creates a bubble of a scalar field given params for initial
//! matter config
class ScalarBubble
{
  public:
    //! A structure for the input params for scalar field properties and initial
    //! conditions
    struct params_t
    {
        double amplitudeSF; //!< Amplitude of bump in initial SF bubble
        std::array<double, AMREX_SPACEDIM>
            centerSF;   //!< Centre of perturbation in initial SF bubble
        double widthSF; //!< Width of bump in initial SF bubble
        double r_zero;  //!< Position of bump relative to centre
    };

    //! The constructor
    ScalarBubble(params_t a_params, double a_dx);

    //! Function to compute the value of all the initial vars on the grid
    template <class data_t>
    void compute(int i, int j, int k, const amrex::Array4<data_t> &state) const;

  protected:
    double m_dx;
    const params_t m_params; //!< The matter initial condition params

    //! Function to compute the value of phi at each point
    template <class data_t>
    data_t compute_phi(Coordinates<data_t> coords) const;
};

#include "ScalarBubble.impl.hpp"

#endif /* SCALARBUBBLE_HPP_ */

/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef BINARYBH_HPP_
#define BINARYBH_HPP_

#include "BoostedBH.hpp"
#include "Cell.hpp"
#include "Coordinates.hpp"
#include "StateVariables.hpp" //This files needs NUM_VARS - total number of components
#include "Tensor.hpp"
#include "simd.hpp"
#include <array>

enum Lapse
{
    ONE,
    PRE_COLLAPSED,
    CHI
};

class BinaryBH
{
  protected:
    double m_dx;
    BoostedBH bh1;
    BoostedBH bh2;
    int m_initial_lapse;

  public:
    // NOLINTBEGIN(bugprone-easily-swappable-parameters)
    BinaryBH(BoostedBH::params_t a_bh1_params, BoostedBH::params_t a_bh2_params,
             double a_dx, int a_initial_lapse = Lapse::PRE_COLLAPSED)
        : m_dx(a_dx), bh1(a_bh1_params), bh2(a_bh2_params),
          m_initial_lapse(a_initial_lapse)
    // NOLINTEND(bugprone-easily-swappable-parameters)
    {
    }

    template <class data_t>
    AMREX_GPU_DEVICE void init_data(int i, int j, int k,
                                    const amrex::CellData<data_t> &cell) const;

  protected:
    template <class data_t>
    AMREX_GPU_DEVICE data_t compute_chi(Coordinates<data_t> coords) const;

    template <class data_t>
    AMREX_GPU_DEVICE Tensor<2, data_t>
    compute_A(data_t chi, Coordinates<data_t> coords) const;
};

#include "BinaryBH.impl.hpp"

#endif /* BINARYBH_HPP_ */

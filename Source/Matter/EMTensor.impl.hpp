/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#if !defined(EMTENSOR_HPP)
#error "This file should only be included through EMTensor.hpp"
#endif

#ifndef EMTENSOR_IMPL_HPP
#define EMTENSOR_IMPL_HPP

#include "CCZ4Geometry.hpp"
#include "Cell.hpp"
#include "DimensionDefinitions.hpp"
#include "FourthOrderDerivatives.hpp"
#include "Interval.hpp"
#include "simd.hpp"

template <class matter_t>
EMTensor<matter_t>::EMTensor(const matter_t &a_matter, const double dx,
                             const int a_c_rho, const Interval a_c_Si,
                             const Interval a_c_Sij)
    : m_matter(a_matter), m_deriv(dx), m_c_rho(a_c_rho), m_c_Si(a_c_Si),
      m_c_Sij(a_c_Sij)
{
    if (m_c_Si.size() != 0)
    {
        // Si is a vector
        AMREX_ASSERT(m_c_Si.size() == DEFAULT_TENSOR_DIM);
    }

    if (m_c_Sij.size() != 0)
    {
        // Sij is a symmetric tensor
        AMREX_ASSERT(m_c_Sij.size() ==
                     DEFAULT_TENSOR_DIM * (DEFAULT_TENSOR_DIM + 1) / 2);
    }
}

template <class matter_t>
template <class data_t>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
EMTensor<matter_t>::compute(int i, int j, int k,
                            const amrex::Array4<data_t> &out_mf,
                            const amrex::Array4<const data_t> &in_mf) const
{
    const auto vars = load_vars<Vars>(in_mf.cellData(i, j, k));
    const auto d1   = m_deriv.template diff1<Vars>(i, j, k, in_mf);

    using namespace TensorAlgebra;

    const auto h_UU  = compute_inverse_sym(vars.h);
    const auto chris = compute_christoffel(d1.h, h_UU);

    const auto emtensor = m_matter.compute_emtensor(vars, d1, h_UU, chris.ULL);

    if (m_c_rho >= 0)
    {
        out_mf(i, j, k, m_c_rho) = emtensor.rho;
    }

    if (m_c_Si.size() > 0)
    {
#if DEFAULT_TENSOR_DIM == 3
        FOR (i)
        {
            out_mf(i, j, k, m_c_Si.begin() + i) = emtensor.Si[i];
        }
#endif
    }

    if (m_c_Sij.size() > 0)
    {
#if DEFAULT_TENSOR_DIM == 3
        out_mf(i, j, k, m_c_Sij.begin())     = emtensor.Sij[0][0];
        out_mf(i, j, k, m_c_Sij.begin() + 1) = emtensor.Sij[0][1];
        out_mf(i, j, k, m_c_Sij.begin() + 2) = emtensor.Sij[0][2];
        out_mf(i, j, k, m_c_Sij.begin() + 3) = emtensor.Sij[1][1];
        out_mf(i, j, k, m_c_Sij.begin() + 4) = emtensor.Sij[1][2];
        out_mf(i, j, k, m_c_Sij.begin() + 5) = emtensor.Sij[2][2];

#endif
    }
}

#endif /* EMTENSOR_IMPL_HPP */

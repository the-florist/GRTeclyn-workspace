/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef SECONDORDERDERIVATIVES_HPP_
#define SECONDORDERDERIVATIVES_HPP_

#include "Cell.hpp"
#include "DimensionDefinitions.hpp"
#include "Tensor.hpp"
#include "simd.hpp"
#include <array>

class SecondOrderDerivatives
{
    private:
        amrex::Real m_dx;
        amrex::Real m_one_over_dx;
        amrex::Real m_one_over_dx2;

    public:
        AMREX_GPU_HOST_DEVICE SecondOrderDerivatives(amrex::Real dx)
        : m_dx(dx), m_one_over_dx(1. / dx), m_one_over_dx2(1. / dx / dx)
        {
        }

        template <class data_t>
        AMREX_GPU_DEVICE AMREX_FORCE_INLINE data_t diff1(const amrex::Real *in_ptr,
                                                         const int idx,
                                                         const int stride) const
        {
            const auto *in = SIMDIFY<data_t>(in_ptr);
            const data_t weight = 0.5;
            
            return (weight * in[idx + stride] -
                    weight * in[idx - stride]) * 
                    m_one_over_dx;
        }

        template <template <typename> class vars_t, class data_t> 
        [[nodiscard]] AMREX_GPU_DEVICE AMREX_FORCE_INLINE auto 
        diff1(int i, int j, int k, const amrex::Array4<const data_t> &state) const 
        {
            vars_t<Tensor<1, data_t>> d1;
            const auto *state_ptr_ijk = state.ptr(i, j, k);
            int j_stride              = static_cast<int>(state.stride.a[0]);
            int k_stride              = static_cast<int>(state.stride.a[1]);
            int n_stride              = static_cast<int>(state.stride.a[2]);

            d1.enum_mapping(
                [&](const int &ivar, Tensor<1, data_t> &var)
                {
                    AMREX_D_TERM(
                        var[0] = diff1<data_t>(state_ptr_ijk + ivar * n_stride, 
                                               0, 1);
                        , 
                        var[1] = diff1<data_t>(state_ptr_ijk + ivar * n_stride, 
                                               0, static_cast<int>(j_stride));
                        ,
                        var[2] = diff1<data_t>(state_ptr_ijk + ivar * n_stride,
                                               0, static_cast<int>(k_stride))
                    );
                });

            return d1;
        }

        template <class data_t>
        AMREX_GPU_DEVICE AMREX_FORCE_INLINE data_t diff2(const amrex::Real *in_ptr,
                                                         const int idx,
                                                         const int stride) const
        {
            const auto *in = SIMDIFY<data_t>(in_ptr);
            
            return (in[idx + stride] +
                    in[idx - stride]
                    - 2. * in[idx]) * 
                    m_one_over_dx2;
        }

        template <template <typename> class vars_t, class data_t> 
        [[nodiscard]] AMREX_GPU_DEVICE AMREX_FORCE_INLINE auto 
        diff2(int i, int j, int k, const amrex::Array4<data_t const> &state) const 
        {
            vars_t<Tensor<2, data_t>> d2;
            const auto *state_ptr_ijk = state.ptr(i, j, k);
            int j_stride              = static_cast<int>(state.stride.a[0]);
            int k_stride              = static_cast<int>(state.stride.a[1]);
            int n_stride              = static_cast<int>(state.stride.a[2]);

            amrex::GpuArray<int, AMREX_SPACEDIM> strides{
                1, static_cast<int>(j_stride),
                static_cast<int>(k_stride)};
            d2.enum_mapping(
                [&](const int &ivar, Tensor<2, data_t> &var)
                {
                    const auto *pvar = state_ptr_ijk + ivar * n_stride;
                    FOR (dir1)
                    {
                        var[dir1][dir1] = diff2<data_t>(pvar, 0, strides[dir1]);
                    }
                });
            return d2;
        }

        template <template <typename> class vars_t, class data_t>
        [[nodiscard]] AMREX_GPU_DEVICE AMREX_FORCE_INLINE auto
        advection(int i, int j, int k, const amrex::Array4<data_t const> &state,
                const Tensor<1, data_t> &vector) const
        {
            vars_t<data_t> advec;
            const auto *state_ptr_ijk = state.ptr(i, j, k);
            int j_stride              = static_cast<int>(state.stride.a[0]);
            int k_stride              = static_cast<int>(state.stride.a[1]);
            int n_stride              = static_cast<int>(state.stride.a[2]);

            amrex::GpuArray<int, AMREX_SPACEDIM> strides{
                1, static_cast<int>(j_stride),
                static_cast<int>(k_stride)};
            advec.enum_mapping(
                [&](const int &ivar, data_t &var)
                {
                    var = 0.;
                    const auto *pvar = state_ptr_ijk + ivar * n_stride;
                    FOR (dir)
                    {
                        var += vector[dir] * diff1<data_t>(pvar, 0, strides[dir]);
                    }
                }
            );
            return advec;
        }

        // Keeping these 3rd order accurate since we don't use them anyways
        template <class data_t>
        AMREX_GPU_DEVICE AMREX_FORCE_INLINE data_t dissipation_term(
            const double *in_ptr, const int idx, const int stride) const
        {
            const auto *const in = SIMDIFY<data_t>(in_ptr);
            data_t weight_vfar   = 1.56250e-2;
            data_t weight_far    = 9.37500e-2;
            data_t weight_near   = 2.34375e-1;
            data_t weight_local  = 3.12500e-1;

            return (weight_vfar * in[idx - 3 * stride] -
                    weight_far * in[idx - 2 * stride] +
                    weight_near * in[idx - stride] - weight_local * in[idx] +
                    weight_near * in[idx + stride] -
                    weight_far * in[idx + 2 * stride] +
                    weight_vfar * in[idx + 3 * stride]) *
                m_one_over_dx;
        }

        template <class data_t, template <typename> class vars_t>
        AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
        add_dissipation(int i, int j, int k, vars_t<data_t> &vars,
                        const amrex::Array4<data_t const> &state,
                        const double factor) const
        {
            const auto *state_ptr_ijk = state.ptr(i, j, k);
            int j_stride              = static_cast<int>(state.stride.a[0]);
            int k_stride              = static_cast<int>(state.stride.a[1]);
            int n_stride              = static_cast<int>(state.stride.a[2]);

            amrex::GpuArray<int, AMREX_SPACEDIM> strides{
                1, static_cast<int>(j_stride),
                static_cast<int>(k_stride)};
            vars.enum_mapping(
                [&](const int &ivar, data_t &var)
                {
                    FOR (dir)
                    {
                        const auto stride  = strides[dir];
                        var               += factor *
                            dissipation_term<data_t>(
                                state_ptr_ijk + ivar * n_stride, 0, stride);
                    }
                });
        }
};

#endif /* SECONDORDERDERIVATIVES__HPP__ */
/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef INFLATONUTILS_HPP_
#define INFLATONUTILS_HPP_

#include "Tensor.hpp"

#include <AMReX_Array.H>
#include <AMReX_Math.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Print.H>

#include <numbers>

#include "InflatonParameters.hpp"

// Trivially-copyable configuration and methods for the stochastic inflaton
// initialisation. The parameters are read once on the host by fill_params()
// and the whole struct is then captured by value into the device kernels, so
// it must stay POD (no amrex::Vector members).
struct InflatonUtils
{
    InflatonParameters m_params;
    InflatonUtils() { m_params.fill_params(); }

    /* Constexpr objects */

    static constexpr amrex::Real tolerance = 1.e-12;

    /* Device-callable functions */

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE static std::uint64_t
    splitmix64(std::uint64_t x)
    {
        x += 0x9E3779B97F4A7C15ULL;
        x  = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ULL;
        x  = (x ^ (x >> 27)) * 0x94D049BB133111EBULL;
        return x ^ (x >> 31);
    }

    // Map 64 random bits to the OPEN interval (0, 1), so that the Rayleigh
    // draw sqrt(-2 ln u) stays finite (u is never exactly 0 or 1).
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE static amrex::Real
    to_unit_open(const std::uint64_t bits)
    {
        return (static_cast<amrex::Real>(bits >> 11) + 0.5) *
               (1.0 / 9007199254740992.0);
    }

    // Nyquist condition
    AMREX_GPU_HOST_DEVICE [[nodiscard]] int flip_index(const int indx) const
    {
        AMREX_ASSERT(m_params.N_fine > 0);
        return amrex::Math::abs(m_params.N_fine - indx);
    }

    // Nyquist condition and calculation of kmag
    AMREX_GPU_HOST_DEVICE [[nodiscard]] int invert_index(const int indx) const
    {
        AMREX_ASSERT(m_params.N_fine > 0);
        return (m_params.N_fine / 2 -
                amrex::Math::abs(m_params.N_fine / 2 - indx));
    }

    // For calculation of polarisation tensors
    AMREX_GPU_HOST_DEVICE [[nodiscard]] int
    invert_index_with_sign(const int indx) const
    {
        AMREX_ASSERT(m_params.N_fine > 0);
        if (indx <= m_params.N_fine / 2)
        {
            return indx;
        }
        return amrex::Math::abs(m_params.N_fine / 2 - indx) -
               m_params.N_fine / 2;
    }

    // Find the magnitude of the Fourier wavevector at this point
    AMREX_GPU_HOST_DEVICE [[nodiscard]] amrex::Real
    get_kmag(amrex::IntVect ivec) const
    {
        AMREX_ASSERT(m_params.L > 0);
        const int k_1 = ivec[0];
        const int k_2 = invert_index(ivec[1]);
        const int k_3 = invert_index(ivec[2]);
        return std::sqrt(k_1 * k_1 + k_2 * k_2 + k_3 * k_3) * 2. *
               amrex::Math::pi<amrex::Real>() / m_params.L;
    }

    // Physical FFT normalisation, shared by the init and extraction classes.
    // CHANGE WITH CARE
    [[nodiscard]] amrex::Real calculate_norm() const
    {
        AMREX_ASSERT(m_params.L > 0);
        return std::pow(
            std::sqrt(2. * amrex::Math::pi<amrex::Real>()) / m_params.L, 3.);
    }

    AMREX_GPU_HOST_DEVICE [[nodiscard]] amrex::Real
    calculate_window_function(const amrex::Real kmag) const
    {
        AMREX_ASSERT(m_params.L > 0 && m_params.Delta > 0);
        const int N_w =
            (m_params.N_coarse != 0) ? m_params.N_coarse : m_params.N;
        const amrex::Real k_cutoff     = std::numbers::sqrt3 * N_w *
                                         amrex::Math::pi<amrex::Real>() /
                                         m_params.L / 5. / 2.;
        const amrex::Real window_width = m_params.L / m_params.Delta;
        return 0.5 * (1.0 - tanh(window_width * (kmag - k_cutoff)));
    }

    //! An orthonormal pair of polarisation basis vectors for a Fourier mode
    struct BasisVectors
    {
        amrex::GpuArray<amrex::Real, 3> mhat;
        amrex::GpuArray<amrex::Real, 3> nhat;
    };

    //! The plus and cross polarisation tensors for a Fourier mode
    struct PolarisationTensors
    {
        Tensor::Rank2 eplus;
        Tensor::Rank2 ecross;
    };

    //! Computes both polarisation tensors for this mode in one call
    AMREX_GPU_HOST_DEVICE [[nodiscard]] PolarisationTensors
    calculate_polarisation_tensors(const amrex::IntVect ivec) const
    {
        // Find basis vectors
        const auto [mhat, nhat] = calculate_basis_vectors(ivec);

        PolarisationTensors pol;
        for (int l = 0; l < 3; l++)
        {
            for (int p = 0; p < 3; p++)
            {
                // Assemble the polarisation tensors
                pol.eplus(l, p)  = mhat[l] * mhat[p] - nhat[l] * nhat[p];
                pol.ecross(l, p) = mhat[l] * nhat[p] + nhat[l] * mhat[p];
            }
        }

        return pol;
    }

    /* Host-only functions */

    // Calculates both basis vectors required for the polarisation tensors
    AMREX_GPU_HOST_DEVICE [[nodiscard]] BasisVectors
    calculate_basis_vectors(const amrex::IntVect ivec) const
    {
        using Vec = amrex::GpuArray<amrex::Real, 3>;

        // Hermitian symmetry inversion on k_2 and k, with sign on the last two
        // indices. (!!) The FT implemented in AMReX symmetrises across the i
        // index, so k_1>= 0 always.
        const auto k_1 = static_cast<amrex::Real>(ivec[0]);
        const auto k_2 =
            static_cast<amrex::Real>(invert_index_with_sign(ivec[1]));
        const auto k_3 =
            static_cast<amrex::Real>(invert_index_with_sign(ivec[2]));

        // Default is the zero mode: mhat = nhat = 0, tensors have no average
        Vec mhat{0., 0., 0.};
        Vec nhat{0., 0., 0.};

        if (ivec == amrex::IntVect{0, 0, 0})
        { /* zero mode: leave as zero */
        }

        else if (k_1 != 0.)
        {
            if (k_2 == 0. && k_3 == 0.)
            {
                mhat = Vec{0., 1., 0.};
                nhat = Vec{0., 0., 1.};
            }
            else
            {
                const amrex::Real i2j2 = k_1 * k_1 + k_2 * k_2;
                const amrex::Real i_j  = std::sqrt(i2j2);
                const amrex::Real norm = std::sqrt(i2j2 * (i2j2 + k_3 * k_3));
                mhat                   = Vec{k_2 / i_j, -k_1 / i_j, 0.};
                nhat =
                    Vec{(k_1 * k_3) / norm, (k_2 * k_3) / norm, -i2j2 / norm};
            }
        }

        else if (k_2 != 0.) // k_1== 0
        {
            if (k_3 == 0.)
            {
                mhat = Vec{0., 0., 1.};
                nhat = Vec{1., 0., 0.};
            }
            else
            {
                const amrex::Real j_k = std::sqrt(k_2 * k_2 + k_3 * k_3);
                mhat                  = Vec{-1., 0., 0.};
                nhat                  = Vec{0., -k_3 / j_k, k_2 / j_k};
            }
        }

        else if (k_3 != 0.) // k_1== 0, k_2 == 0
        {
            mhat = Vec{1., 0., 0.};
            nhat = Vec{0., 1., 0.};
        }

        else
        {
            // Unreachable: (i, j, k) == (0, 0, 0) only when iv == {0,0,0},
            // which is handled by the zero-mode branch above.
            AMREX_ASSERT_WITH_MESSAGE(false,
                                      "InflatonUtils::calculate_basis_vectors, "
                                      "Fourier grid point not covered.");
        }

        // Apply the internal rotation in the +/x decomposition basis, if
        // requested
        if (m_params.alpha != 0.)
        {
            const amrex::Real alpha_rad =
                m_params.alpha * amrex::Math::pi<amrex::Real>() / 180.;
            const amrex::Real cos_alpha = std::cos(alpha_rad);
            const amrex::Real sin_alpha = std::sin(alpha_rad);
            Vec mhat_rot;
            Vec nhat_rot;
            for (int l = 0; l < 3; l++)
            {
                mhat_rot[l] = cos_alpha * mhat[l] + sin_alpha * nhat[l];
                nhat_rot[l] = -sin_alpha * mhat[l] + cos_alpha * nhat[l];
            }
            mhat = mhat_rot;
            nhat = nhat_rot;
        }

        return {.mhat = mhat, .nhat = nhat};
    }

    // Applies above Nyquist conditions to a given MF
    void apply_nyquist_conditions(amrex::cMultiFab &field)
    {
        AMREX_ASSERT(m_params.N_fine > 0);

        // Slice to the POD base so the kernel captures config by value, not via
        // the (host) this pointer
        const InflatonUtils cfg = *this;

        int num_components = field.nComp();
        for (amrex::MFIter mfi(field); mfi.isValid(); ++mfi)
        {
            // The geometry for this MPI rank
            const amrex::Box &fabbox = mfi.fabbox();
            amrex::Array4<amrex::GpuComplex<amrex::Real>> const &field_ptr =
                field.array(mfi);

            amrex::ParallelFor(
                fabbox,
                [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
                {
                    amrex::IntVect iv = {i, j, k};

                    if ((i == 0 || i == cfg.m_params.N_fine / 2) &&
                        (j == 0 || j == cfg.m_params.N_fine / 2) &&
                        (k == 0 || k == cfg.m_params.N_fine / 2))
                    {
                        for (int comp = 0; comp < num_components; comp++)
                        {
                            amrex::GpuComplex<amrex::Real> temp(
                                field_ptr(i, j, k, comp).real(), 0.);
                            field_ptr(i, j, k, comp) = temp;
                        }
                    }

                    else if (i == 0 || i == cfg.m_params.N_fine / 2)
                    {
                        const int n_half = cfg.m_params.N_fine / 2;
                        if ((k > n_half && j == n_half) ||
                            (k == 0 && j > n_half) || (k > n_half && j == 0) ||
                            (k == n_half && j > n_half))
                        {
                            for (int comp = 0; comp < num_components; comp++)
                            {
                                amrex::GpuComplex<amrex::Real> temp(
                                    field_ptr(i, cfg.invert_index(j),
                                              cfg.invert_index(k), comp)
                                        .real(),
                                    -field_ptr(i, cfg.invert_index(j),
                                               cfg.invert_index(k), comp)
                                         .imag());
                                field_ptr(i, j, k, comp) = temp;
                            }
                        }

                        else if (j > n_half)
                        {
                            for (int comp = 0; comp < num_components; comp++)
                            {
                                amrex::GpuComplex<amrex::Real> temp(
                                    field_ptr(i, cfg.invert_index(j),
                                              cfg.flip_index(k), comp)
                                        .real(),
                                    -field_ptr(i, cfg.invert_index(j),
                                               cfg.flip_index(k), comp)
                                         .imag());
                                field_ptr(i, j, k, comp) = temp;
                            }
                        }
                    }
                });

            amrex::Gpu::streamSynchronize();
        }
    }
};

#endif /* INFLATONUTILS_HPP_ */

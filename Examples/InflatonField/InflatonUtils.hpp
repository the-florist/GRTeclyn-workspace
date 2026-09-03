/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef INFLATONUTILS_HPP_
#define INFLATONUTILS_HPP_

#include "GRParmParse.hpp"
#include "Tensor.hpp"

#include <AMReX_Array.H>
#include <AMReX_Math.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Print.H>
#include <AMReX_Vector.H>

#include "Potential.hpp"

struct InflatonParameters
{
    // Parameters, populated once on the host by fill_params()
    amrex::Real Mp{1.};
    amrex::Real L{0.};
    amrex::Real Delta{1.};
    amrex::Real alpha{0.};
    int N{0};
    int N_fine{0};
    int N_coarse{0};
    int scalar_init{0};
    int tensor_init{0};
    int use_rand{1};
    int use_window{0};
    int test_normalisation{0};
    int random_seed{3539263};

    amrex::Real phi0{0.};      //!< Background scalar-field value
    amrex::Real H0{0.};        //!< Initial Hubble parameter
    amrex::Real epsilon_1{0.}; //!< First slow-roll parameter
    amrex::Real epsilon_2{0.}; //!< Second slow-roll parameter

    // Read the scalar parameters once from the global GRParmParse table
    void fill_params()
    {
        GRParmParse pp;
        amrex::Vector<int> ncell;
        pp.getarr("amr.n_cell", ncell);
        N      = ncell[0];
        N_fine = N;
        pp.query("N_fine", N_fine);
        N_coarse = N;
        pp.query("N_coarse", N_coarse);
        amrex::Vector<amrex::Real> prob_extent;
        pp.getarr("geometry.prob_extent", prob_extent);
        L = prob_extent[0];

        amrex::Real G_Newton = 1.;
        pp.query("scalar_field.G_Newton", G_Newton);
        Mp              = 1. / std::sqrt(G_Newton);
        amrex::Real Pi0 = 0., V = 0., dV = 0.;

        GRParmParse init_pp("init");
        init_pp.query("scalar_init", scalar_init);
        init_pp.query("tensor_init", tensor_init);
        init_pp.query("use_rand", use_rand);
        init_pp.query("use_window", use_window);
        init_pp.query("test_normalisation", test_normalisation);
        init_pp.query("alpha", alpha);
        init_pp.query("Delta", Delta);
        init_pp.query("random_seed", random_seed);
        init_pp.get("background_phi", phi0);
        init_pp.query("background_dphi", Pi0);

        Potential potential;
        potential.compute_background_potential(V, dV, phi0);

        H0 = calculate_H0(G_Newton, Pi0, V);

        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
            scalar_init == 0 || Pi0 != 0.,
            "Config::fill_params, init.background_dphi must be "
            "non-zero when init.scalar_init is enabled");
        if (Pi0 != 0.)
        {
            epsilon_1 = std::pow(Pi0 / H0, 2.) / 2. / std::pow(Mp, 2.);
            epsilon_2 =
                6. + dV / Pi0 / H0 - std::pow(Pi0 / H0, 2.) / std::pow(Mp, 2.);
        }
    }

    void check_params(const amrex::Vector<int> &ncell,
                      const amrex::Vector<amrex::Real> &prob_extent) const
    {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
            Mp != 0., "Config::check_params, Mp must be non-zero");
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
            N_fine >= N, "Config::check_params, N_fine must be >= N");
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
            N_coarse <= N, "Config::check_params, N_coarse must be <= N");
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
            L != 0., "Config::check_params, L must be non-zero");
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
            Delta != 0., "Config::check_params, Delta must be non-zero");

        for (int d = 1; d < static_cast<int>(ncell.size()); d++)
        {
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
                ncell[d] == ncell[0],
                "Config::check_params, grid must be cubic "
                "(amr.n_cell must be equal in all directions)");
        }
        for (int d = 1; d < static_cast<int>(prob_extent.size()); d++)
        {
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
                prob_extent[d] == prob_extent[0],
                "Config::check_params, domain must be cubic "
                "(geometry.prob_extent must be equal in all directions)");
        }
    }

    AMREX_GPU_HOST_DEVICE static amrex::Real
    calculate_H0(amrex::Real G, amrex::Real Pi, amrex::Real V)
    {
        return sqrt((8. * amrex::Math::pi<amrex::Real>() * G / 3.) *
                    (0.5 * pow(Pi, 2.) + V));
    }
};

// Trivially-copyable configuration and methods for the stochastic inflaton
// initialisation. The parameters are read once on the host by fill_params()
// and the whole struct is then captured by value into the device kernels, so
// it must stay POD (no amrex::Vector members).
struct InflatonUtils
{
    InflatonParameters m_params;
    InflatonUtils() { m_params.fill_params(); }

    /* Constexpr objects */

    // Look-up table
    // Used to construct polarisation basis tensors
    static constexpr int look_up_table[3][3]{
        {0, 1, 2},
        {1, 3, 4},
        {2, 4, 5}
    };
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
    AMREX_GPU_HOST_DEVICE inline int flip_index(const int indx) const
    {
        AMREX_ASSERT(m_params.N_fine > 0);
        return amrex::Math::abs(m_params.N_fine - indx);
    }

    // Nyquist condition and calculation of kmag
    AMREX_GPU_HOST_DEVICE inline int invert_index(const int indx) const
    {
        AMREX_ASSERT(m_params.N_fine > 0);
        return (int)(m_params.N_fine / 2 -
                     amrex::Math::abs(m_params.N_fine / 2 - indx));
    }

    // For calculation of polarisation tensors
    AMREX_GPU_HOST_DEVICE inline int
    invert_index_with_sign(const int indx) const
    {
        AMREX_ASSERT(m_params.N_fine > 0);
        if (indx <= m_params.N_fine / 2)
        {
            return indx;
        }
        else
        {
            return amrex::Math::abs(m_params.N_fine / 2 - indx) -
                   m_params.N_fine / 2;
        }
    }

    // Find the magnitude of the Fourier wavevector at this point
    AMREX_GPU_HOST_DEVICE inline amrex::Real get_kmag(amrex::IntVect iv) const
    {
        AMREX_ASSERT(m_params.L > 0);
        const int i = iv[0];
        const int j = invert_index(iv[1]);
        const int k = invert_index(iv[2]);
        return std::sqrt(i * i + j * j + k * k) * 2. *
               amrex::Math::pi<amrex::Real>() / m_params.L;
    }

    // Physical FFT normalisation, shared by the init and extraction classes.
    // CHANGE WITH CARE
    inline amrex::Real calculate_norm() const
    {
        AMREX_ASSERT(m_params.L > 0);
        return std::pow(
            std::sqrt(2. * amrex::Math::pi<amrex::Real>()) / m_params.L, 3.);
    }

    AMREX_GPU_HOST_DEVICE inline amrex::Real
    calculate_window_function(const amrex::Real kmag) const
    {
        AMREX_ASSERT(m_params.L > 0 && m_params.Delta > 0);
        const int N_w =
            (m_params.N_coarse != 0) ? m_params.N_coarse : m_params.N;
        const amrex::Real ks = std::sqrt(3.) * N_w *
                               amrex::Math::pi<amrex::Real>() / m_params.L /
                               5. / 2.;
        const amrex::Real Dt = m_params.L / m_params.Delta;
        return 0.5 * (1.0 - tanh(Dt * (kmag - ks)));
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
    AMREX_GPU_HOST_DEVICE inline PolarisationTensors
    calculate_polarisation_tensors(const amrex::IntVect iv) const
    {
        // Find basis vectors
        const auto [mhat, nhat] = calculate_basis_vectors(iv);

        PolarisationTensors pol;
        for (int l = 0; l < 3; l++)
            for (int p = 0; p < 3; p++)
            {
                // Assemble the polarisation tensors
                pol.eplus(l, p)  = mhat[l] * mhat[p] - nhat[l] * nhat[p];
                pol.ecross(l, p) = mhat[l] * nhat[p] + nhat[l] * mhat[p];
            }

        return pol;
    }

    /* Host-only functions */

    // Calculates both basis vectors required for the polarisation tensors
    AMREX_GPU_HOST_DEVICE inline BasisVectors
    calculate_basis_vectors(const amrex::IntVect iv) const
    {
        using Vec = amrex::GpuArray<amrex::Real, 3>;

        // Hermitian symmetry inversion on j and k, with sign on the last two
        // indices. (!!) The FT implemented in AMReX symmetrises across the i
        // index, so i >= 0 always.
        const amrex::Real i = static_cast<amrex::Real>(iv[0]);
        const amrex::Real j =
            static_cast<amrex::Real>(invert_index_with_sign(iv[1]));
        const amrex::Real k =
            static_cast<amrex::Real>(invert_index_with_sign(iv[2]));

        // Default is the zero mode: mhat = nhat = 0, tensors have no average
        Vec mhat{0., 0., 0.};
        Vec nhat{0., 0., 0.};

        if (iv == amrex::IntVect{0, 0, 0})
        { /* zero mode: leave as zero */
        }

        else if (i != 0.)
        {
            if (j == 0. && k == 0.)
            {
                mhat = Vec{0., 1., 0.};
                nhat = Vec{0., 0., 1.};
            }
            else
            {
                const amrex::Real i2j2 = i * i + j * j;
                const amrex::Real ij   = std::sqrt(i2j2);
                const amrex::Real n    = std::sqrt(i2j2 * (i2j2 + k * k));
                mhat                   = Vec{j / ij, -i / ij, 0.};
                nhat = Vec{(i * k) / n, (j * k) / n, -i2j2 / n};
            }
        }

        else if (j != 0.) // i == 0
        {
            if (k == 0.)
            {
                mhat = Vec{0., 0., 1.};
                nhat = Vec{1., 0., 0.};
            }
            else
            {
                const amrex::Real jk = std::sqrt(j * j + k * k);
                mhat                 = Vec{-1., 0., 0.};
                nhat                 = Vec{0., -k / jk, j / jk};
            }
        }

        else if (k != 0.) // i == 0, j == 0
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
            const amrex::Real a =
                m_params.alpha * amrex::Math::pi<amrex::Real>() / 180.;
            const amrex::Real ca = std::cos(a);
            const amrex::Real sa = std::sin(a);
            Vec mp, np;
            for (int l = 0; l < 3; l++)
            {
                mp[l] = ca * mhat[l] + sa * nhat[l];
                np[l] = -sa * mhat[l] + ca * nhat[l];
            }
            mhat = mp;
            nhat = np;
        }

        return {mhat, nhat};
    }

    // Applies above Nyquist conditions to a given MF
    inline void apply_nyquist_conditions(amrex::cMultiFab &field)
    {
        AMREX_ASSERT(m_params.N_fine > 0);

        // Slice to the POD base so the kernel captures config by value, not via
        // the (host) this pointer
        const InflatonUtils cfg = *this;

        int num_components = field.nComp();
        for (amrex::MFIter mfi(field); mfi.isValid(); ++mfi)
        {
            // The geometry for this MPI rank
            const amrex::Box &bx = mfi.fabbox();
            amrex::Array4<amrex::GpuComplex<amrex::Real>> const &field_ptr =
                field.array(mfi);

            amrex::ParallelFor(
                bx,
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
                        const int nh = cfg.m_params.N_fine / 2;
                        if ((k > nh && j == nh) || (k == 0 && j > nh) ||
                            (k > nh && j == 0) || (k == nh && j > nh))
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

                        else if (j > nh)
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

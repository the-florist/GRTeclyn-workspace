/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef CONFIG_HPP_
#define CONFIG_HPP_

#include "GRParmParse.hpp"
#include "Tensor.hpp"

#include <AMReX_Array.H>
#include <AMReX_Math.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Print.H>
#include <AMReX_Vector.H>

#include "Potential.hpp"
#include "TensorTests.hpp"
#include "Utils.hpp"

// Trivially-copyable configuration and methods for the stochastic inflaton
// initialisation. The parameters are read once on the host by fill_params()
// and the whole struct is then captured by value into the device kernels, so
// it must stay POD (no amrex::Vector members).
struct Config
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

    /* Device-callable functions */

    // Nyquist condition
    AMREX_GPU_HOST_DEVICE inline int flip_index(const int indx) const
    {
        AMREX_ASSERT(N_fine > 0);
        return amrex::Math::abs(N_fine - indx);
    }

    // Nyquist condition and calculation of kmag
    AMREX_GPU_HOST_DEVICE inline int invert_index(const int indx) const
    {
        AMREX_ASSERT(N_fine > 0);
        return (int)(N_fine / 2 - amrex::Math::abs(N_fine / 2 - indx));
    }

    // For calculation of polarisation tensors
    AMREX_GPU_HOST_DEVICE inline int
    invert_index_with_sign(const int indx) const
    {
        AMREX_ASSERT(N_fine > 0);
        if (indx <= N_fine / 2)
        {
            return indx;
        }
        else
        {
            return amrex::Math::abs(N_fine / 2 - indx) - N_fine / 2;
        }
    }

    // Find the magnitude of the Fourier wavevector at this point
    AMREX_GPU_HOST_DEVICE inline amrex::Real get_kmag(amrex::IntVect iv) const
    {
        AMREX_ASSERT(L > 0);
        const int i = iv[0];
        const int j = invert_index(iv[1]);
        const int k = invert_index(iv[2]);
        return std::sqrt(i * i + j * j + k * k) * 2. *
               amrex::Math::pi<amrex::Real>() / L;
    }

    // Physical FFT normalisation, shared by the init and extraction classes.
    // CHANGE WITH CARE
    inline amrex::Real calculate_norm() const
    {
        AMREX_ASSERT(L > 0);
        return std::pow(std::sqrt(2. * amrex::Math::pi<amrex::Real>()) / L, 3.);
    }

    AMREX_GPU_HOST_DEVICE inline amrex::Real
    calculate_window_function(const amrex::Real kmag) const
    {
        AMREX_ASSERT(L > 0 && Delta > 0);
        const int N_w = (N_coarse != 0) ? N_coarse : N;
        const amrex::Real ks =
            std::sqrt(3.) * N_w * amrex::Math::pi<amrex::Real>() / L / 5. / 2.;
        const amrex::Real Dt = L / Delta;
        return 0.5 * (1.0 - tanh(Dt * (kmag - ks)));
    }

    AMREX_GPU_HOST_DEVICE static amrex::Real
    calculate_H0(amrex::Real G, amrex::Real Pi, amrex::Real V)
    {
        return sqrt((8. * amrex::Math::pi<amrex::Real>() * G / 3.) *
                    (0.5 * pow(Pi, 2.) + V));
    }

    //! An orthonormal pair of polarisation basis vectors for a Fourier mode
    struct BasisVectors
    {
        amrex::GpuArray<amrex::Real, 3> mhat;
        amrex::GpuArray<amrex::Real, 3> nhat;
    };

    //! Computes both polarisation basis vectors for this mode in one call
    AMREX_GPU_HOST_DEVICE inline BasisVectors
    calculate_basis_vectors(const amrex::IntVect iv) const;

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

    inline void test_polarisation_normalisation(const amrex::cMultiFab &kfield)
    {
        for (amrex::MFIter mfi(kfield); mfi.isValid(); ++mfi)
        {
            const amrex::Box &bx    = mfi.fabbox();
            const amrex::IntVect lo = bx.smallEnd();
            const amrex::IntVect hi = bx.bigEnd();
            for (int k = lo[2]; k <= hi[2]; k++)
                for (int j = lo[1]; j <= hi[1]; j++)
                    for (int i = lo[0]; i <= hi[0]; i++)
                    {
                        const amrex::IntVect iv{i, j, k};
                        const auto [mhat, nhat] = calculate_basis_vectors(iv);
                        TensorTests::test_vector_orthonorm(iv, mhat, nhat);

                        const auto [eplus, ecross] =
                            calculate_polarisation_tensors(iv);
                        TensorTests::test_polarisation_tensor_orthonorm(
                            iv, eplus, ecross);
                    }
        }
    }

    // Applies above Nyquist conditions to a given MF
    inline void apply_nyquist_conditions(amrex::cMultiFab &field);
};

#include "Config.impl.hpp"

#endif /* CONFIG_HPP_ */

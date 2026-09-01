/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef INFLATIONCONFIG_HPP_
#define INFLATIONCONFIG_HPP_

#include "GRParmParse.hpp"
#include "Tensor.hpp"

#include <AMReX_MultiFab.H>
#include <AMReX_Print.H>
#include <AMReX_Vector.H>
#include <AMReX_Array.H>
#include <AMReX_Math.H>

#include "InflationUtils.hpp"
#include "TensorTests.hpp"

// Trivially-copyable configuration and methods for the stochastic inflaton
// initialisation. The parameters are read once on the host by fill_params()
// and the whole struct is then captured by value into the device kernels, so
// it must stay POD (no amrex::Vector members).
struct InflationConfig
{
    // Parameters, populated once on the host by fill_params()
    amrex::Real Mp{1.};
    amrex::Real L{0.};
    amrex::Real Delta{1.};
    amrex::Real alpha{0.};
    int N{0};
    int N_fine{0};
    int N_coarse{0};
    int read_from_stoiic{0};
    int scalar_init{0};
    int tensor_init{0};
    int use_rand{1};
    int use_window{0};
    int test_normalisation{0};
    int random_seed{3539263};

    // STOIIC spectra uploaded to the device, flattened row-major as
    // [row][mode]. Set by RandomFieldInit; null when STOIIC is unused.
    int n_modes{0};
    const amrex::Real *init_k_ptr{nullptr};
    const amrex::Real *scalar_ps_ptr{nullptr};
    const amrex::Real *tensor_ps_ptr{nullptr};

    // Read the scalar parameters once from the global GRParmParse table
    void fill_params()
    {
        GRParmParse pp;
        pp.get("N", N);
        N_fine = N;
        pp.query("N_fine", N_fine);
        N_coarse = N;
        pp.query("N_coarse", N_coarse);
        pp.get("L", L);

        amrex::Real G_Newton = 1.;
        pp.query("G_Newton", G_Newton);
        Mp = 1. / std::sqrt(G_Newton);

        GRParmParse randominit_pp("randominit");
        randominit_pp.query("read_from_STOIIC", read_from_stoiic);
        randominit_pp.query("scalar_init", scalar_init);
        randominit_pp.query("tensor_init", tensor_init);
        randominit_pp.query("use_rand", use_rand);
        randominit_pp.query("use_window", use_window);
        randominit_pp.query("test_normalisation", test_normalisation);
        randominit_pp.query("alpha", alpha);
        randominit_pp.query("Delta", Delta);
        randominit_pp.query("random_seed", random_seed);
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
        return (int)(N_fine/2 - amrex::Math::abs(N_fine/2 - indx));
    }

    // For calculation of polarisation tensors
    AMREX_GPU_HOST_DEVICE inline int invert_index_with_sign(const int indx) const
    {
        AMREX_ASSERT(N_fine > 0);
        if (indx <= N_fine/2) { return indx; }
        else { return amrex::Math::abs(N_fine/2 - indx) - N_fine/2; }
    }

    // Find the magnitude of the Fourier wavevector at this point
    AMREX_GPU_HOST_DEVICE inline amrex::Real get_kmag(amrex::IntVect iv) const
    {
        AMREX_ASSERT(L > 0);
        const int i = iv[0];
        const int j = invert_index(iv[1]);
        const int k = invert_index(iv[2]);
        return std::sqrt(i*i + j*j + k*k) * 2. * amrex::Math::pi<amrex::Real>() / L;
    }

    // Physical FFT normalisation, shared by the init and extraction classes.
    // CHANGE WITH CARE
    inline amrex::Real norm() const
    {
        AMREX_ASSERT(L > 0);
        return std::pow(std::sqrt(2. * amrex::Math::pi<amrex::Real>()) / L, 3.);
    }

    AMREX_GPU_HOST_DEVICE inline amrex::Real
    window_function(const amrex::Real kmag) const
    {
        AMREX_ASSERT(L > 0 && Delta > 0);
        const int N_w = (N_coarse != 0) ? N_coarse : N;
        const amrex::Real ks =
            std::sqrt(3.) * N_w * amrex::Math::pi<amrex::Real>() / L / 5. / 2.;
        const amrex::Real Dt = L / Delta;
        return 0.5 * (1.0 - tanh(Dt * (kmag - ks)));
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
        Tensor<2, amrex::Real> eplus;
        Tensor<2, amrex::Real> ecross;
    };

    //! Computes both polarisation tensors for this mode in one call
    AMREX_GPU_HOST_DEVICE inline PolarisationTensors
    calculate_polarisation_tensors(const amrex::IntVect iv) const
    {
        // Find basis vectors
        const auto [mhat, nhat] = calculate_basis_vectors(iv);

        PolarisationTensors pol;
        for (int l=0; l<3; l++) for (int p=0; p<3; p++)
        {
            // Assemble the polarisation tensors
            pol.eplus[l][p] = mhat[l]*mhat[p] - nhat[l]*nhat[p];
            pol.ecross[l][p] = mhat[l]*nhat[p] + nhat[l]*mhat[p];
        }

        return pol;
    }

    /* Host-only functions */

    inline void test_polarisation_normalisation(const amrex::cMultiFab &kfield)
    {
        for (amrex::MFIter mfi(kfield); mfi.isValid(); ++mfi)
        {
            const amrex::Box &bx = mfi.fabbox();
            const amrex::IntVect lo = bx.smallEnd();
            const amrex::IntVect hi = bx.bigEnd();
            for (int k = lo[2]; k <= hi[2]; k++)
            for (int j = lo[1]; j <= hi[1]; j++)
            for (int i = lo[0]; i <= hi[0]; i++)
            {
                const amrex::IntVect iv{i, j, k};
                const auto [mhat, nhat] = calculate_basis_vectors(iv);
                TensorTests::Test_vector_orthonorm(iv, mhat, nhat);

                const auto [eplus, ecross] = calculate_polarisation_tensors(iv);
                TensorTests::Test_polarisation_tensor_orthonorm(iv, eplus, ecross);
            }
        }
    }

    // Applies above Nyquist conditions to a given MF
    inline void apply_nyquist_conditions(amrex::cMultiFab &field);
};

#include "InflationConfig.impl.hpp"

#endif /* INFLATIONCONFIG */

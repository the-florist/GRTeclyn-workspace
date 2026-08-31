/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef INFLATIONCONFIG_HPP_
#define INFLATIONCONFIG_HPP_

#include "Tensor.hpp"

#include <AMReX_MultiFab.H>
#include <AMReX_Print.H>
#include <AMReX_Vector.H>
#include <AMReX_Array.H>

#include "InflationUtils.hpp"
#include "TensorTests.hpp"

struct InflationParams
{
    /* Shared parameters */
    amrex::Real Mp{1.};             //!< Energy scale of the problem

    // Basic initialisation flags
    int read_from_stoiic{0};   //!< Whether to read spectrum from stoiic input
    int tensor_init{0};        //!< Determines whether tensor perturbations run
    int scalar_init{0};        //!< Read in perturbations from STOIIC dparams
    int use_rand{1};           //!< Choose whether to use random initial conditions
    int use_window{0};         //!< Choose whether to use window function
    int test_normalisation{0}; //!< Run the one-time polarisation-basis
                               //!< orthonormality check

    // Grid parameters
    amrex::Real L{0};       //!< Length of the box
    int N{0};               //!< Grid resolution (number of points per dimension)
    int N_fine{0};          //!< Fine resolution to downsample from,
                            //!< used for convergence testing
    int N_coarse{0};        //!< Coarse resolution to use for the cutoff mode
                            //!< Used for convergence testing.

    // Field construction parameters
    amrex::Real A{1.};         //!< Amplitude factor (for basic tests)
    int random_seed{3539263};  //!< Seed for random number generator
    amrex::Real alpha{0.};     //!< Internal rotation angle in the +/x
                               //!< decomposition basis
    amrex::Real Delta{1.};     //!< window's width, measured like L/Delta

    // Extraction parameters
    int calc_binned_power_spectrum{0};   //!< Choose whether to extract the
                                         //!< binned power spectrum
    int plot_int{0};
    int calc_higher_order_statistics{0}; //!< Choose whether to print higher-order
                                         //!< statistics on the fields
    int num_orders{0};                   //!< Number of moments to print
                                         //!< (required by vector read-in)

    // STOIIC spectra uploaded to the device, flattened row-major as
    // [row][mode]. Set by RandomFieldInit; null when STOIIC is unused.
    int n_modes{0};                             //!< Number of k modes read in
    const amrex::Real *init_k_ptr{nullptr};     //!< k values, length n_modes
    const amrex::Real *scalar_ps_ptr{nullptr};  //!< 8 x n_modes scalar spectra
    const amrex::Real *tensor_ps_ptr{nullptr};  //!< 4 x n_modes tensor spectra

    /* Shared device-callable functions */

    // Nyquist condition
    AMREX_GPU_HOST_DEVICE inline int flip_index(const int indx) const
    {
        AMREX_ASSERT(N_fine > 0);
        return std::abs(N_fine - indx);
    }

    // Nyquist condition and calculation of kmag
    AMREX_GPU_HOST_DEVICE inline int invert_index(const int indx) const
    {
        AMREX_ASSERT(N_fine > 0);
        return (int)(N_fine/2 - std::abs(N_fine/2 - indx));
    }

    // For calculation of polarisation tensors
    AMREX_GPU_HOST_DEVICE inline int invert_index_with_sign(const int indx) const
    {
        AMREX_ASSERT(N_fine > 0);
        if (indx <= N_fine/2) { return indx; }
        else { return std::abs(N_fine/2 - indx) - N_fine/2; }
    }

    // Find the magnitude of the Fourier wavevector at this point
    AMREX_GPU_HOST_DEVICE inline amrex::Real get_kmag(amrex::IntVect iv) const
    {
        AMREX_ASSERT(L > 0);
        const int i = iv[0];
        const int j = invert_index(iv[1]);
        const int k = invert_index(iv[2]);
        return std::sqrt(i*i + j*j + k*k) * 2. * M_PI / L;
    }

    // Physical FFT normalisation, shared by the init and extraction classes.
    // CHANGE WITH CARE
    inline amrex::Real norm() const
    {
        AMREX_ASSERT(L > 0);
        return std::pow(std::sqrt(2. * M_PI) / L, 3.);
    }

    AMREX_GPU_HOST_DEVICE inline amrex::Real
    window_function(const amrex::Real kmag) const
    {
        AMREX_ASSERT(L > 0 && Delta > 0);
        const int N_w = (N_coarse != 0) ? N_coarse : N;
        const amrex::Real ks = std::sqrt(3.) * N_w * M_PI / L / 5. / 2.;
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
    calculate_basis_vectors(const amrex::IntVect iv);

    //! The plus and cross polarisation tensors for a Fourier mode
    struct PolarisationTensors
    {
        Tensor<2, amrex::Real> eplus;
        Tensor<2, amrex::Real> ecross;
    };

    //! Computes both polarisation tensors for this mode in one call
    AMREX_GPU_HOST_DEVICE inline PolarisationTensors
    calculate_polarisation_tensors(const amrex::IntVect iv)
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
};

struct InflationConfig : public InflationParams
{
    amrex::Vector<int> orders;           //!< Moment orders to print for
                                         //!< extracted fields

    // STOIIC read-in structures
    //!< ks printed by STOIIC, at which Fourier-space fields are provided
    amrex::Vector<amrex::Real> init_k;
    //!< Structure: four fields * two components, power spec values
    amrex::Vector<amrex::Vector<amrex::Real>> scalar_ps;
    //!< Structure: two fields * two components, power spec values
    amrex::Vector<amrex::Vector<amrex::Real>> tensor_ps;

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

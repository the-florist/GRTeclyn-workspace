/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef INFLATIONCONFIG_HPP_
#define INFLATIONCONFIG_HPP_

#include "Cell.hpp"
#include "InitialScalarData.hpp"
#include "VarsTools.hpp"
#include "FilesystemTools.hpp"
#include "Potential.hpp"
#include "Tensor.hpp"
#include <fstream>
#include <random>

#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>
#include <AMReX_FFT.H>
#include <AMReX_Random.H>
#include <AMReX_Print.H>
#include <AMReX_Vector.H>
#include <AMReX_Array.H>

#include "InflationUtils.hpp"
#include "TensorTests.hpp"

struct InflationConfig
{
    /* Shared parameters */
    amrex::Real Mp = 1.;             //!< Energy scale of the problem

    // Basic initialisation flags
    int read_from_stoiic = 0;   //!< Whether to read spectrum from stoiic dparams.txt input
    int tensor_init = 0;        //!< Determines whether tensor perturbations are calculated
    int scalar_init = 0;  //!< Read in perturbations from STOIIC dparams
    int use_rand = 1;           //!< Choose whether to use random initial conditions
    int use_window = 0;         //!< Choose whether to use window function

    // Grid parameters
    amrex::Real L = 0;                   //!< Length of the box
    int N = 0;               //!< used to read in the private N variable
    int N_fine = 0;                 //!< Fine resolution to downsample from, 
                                //!< used for convergence testing
    int N_coarse = 0;           //!< Coarse resolution to use for the cutoff mode
                                //!< Used for convergence testing.

    // Field construction parameters
    amrex::Real A = 0;                   //!< Amplitude factor (for basic tests)
    int random_seed = 3539263;  //!< Seed for random number generator
    amrex::Real alpha = 0.;          //!< Internal rotation angle in the +/x decomposition basis
    amrex::Real kstar = 0;               //!< window's cut-off mode, measured in units of 2pi/L
    amrex::Real Delta = 0;               //!< window's width, measured like L/Delta

    // Extraction parameters
    int calc_binned_power_spectrum = 0;   //!< Choose whether to extract the binned power spectrum
    int plot_int = 0;
    int calc_higher_order_statistics = 0; //!< Choose whether to print higher-order statistics on the fields
    int num_orders = 0;                       //!< Number of moments to print (required by vector read-in)
    amrex::Vector<int> orders;                   //!< Moment orders to print for extracted fields

    // STOIIC read-in structures
    amrex::Vector<amrex::Real> init_k;                  //!< ks printed by STOIIC, at which Fourier-space fields are provided
    amrex::Vector<amrex::Vector<amrex::Real>> scalar_ps;       //!< Structure: four fields * two components, power spec values
    amrex::Vector<amrex::Vector<amrex::Real>> tensor_ps;       //!< Structure: two fields * two components, power spec values

    /* Shared functions */

    // Nyquist condition
    inline int flip_index(const int indx) 
    {
        AMREX_ASSERT(N_fine > 0); 
        return std::abs(N_fine - indx); 
    }

    // Nyquist condition and calculation of kmag
    inline int invert_index(const int indx) 
    { 
        AMREX_ASSERT(N_fine > 0);
        return (int)(N_fine/2 - std::abs(N_fine/2 - indx)); 
    }

    // For calculation of polarisation tensors
    inline int invert_index_with_sign(const int indx) 
    { 
        AMREX_ASSERT(N_fine > 0);
        if (indx <= N_fine/2) { return indx; }
        else { return std::abs(N_fine/2 - indx) - N_fine/2; }
    }

    // Find the magnitude of the Fourier wavevector at this point
    inline amrex::Real get_kmag(amrex::IntVect iv)
    {
        AMREX_ASSERT(L > 0);
        const int i = iv[0];
        const int j = invert_index(iv[1]);
        const int k = invert_index(iv[2]);
        return std::sqrt(i*i + j*j + k*k) * 2. * M_PI / L;
    }

    inline amrex::GpuArray<amrex::Real, 3> 
    calculate_basis_vector(const amrex::IntVect iv, 
                           const int which_vector);

    inline Tensor<2, amrex::Real> 
    calculate_polarisation_tensor(const amrex::IntVect iv,
                                  const int which_pol)
    {
        // Find basis vectors
        amrex::GpuArray<amrex::Real, 3> mhat = calculate_basis_vector(iv, 0);
        amrex::GpuArray<amrex::Real, 3> nhat = calculate_basis_vector(iv, 1);
        TensorTests::Test_vector_orthonorm(iv, mhat, nhat);

        Tensor<2, amrex::Real> eplus, ecross; 
        for (int l=0; l<3; l++) for (int p=0; p<3; p++)
        {
            // Assemble the polarisation tensors
            eplus[l][p] = mhat[l]*mhat[p] - nhat[l]*nhat[p];
            ecross[l][p] = mhat[l]*nhat[p] + nhat[l]*mhat[p];
        }

        TensorTests::Test_polarisation_tensor_orthonorm(iv, eplus, ecross);

        if (which_pol == 0) { return eplus; }
        else if (which_pol == 1) { return ecross; }
        else 
        {
            amrex::Error("InflationConfig::calculate_polarisation_tensor, "
                         "polarisation flag is not set correctly.");
            return eplus;
        }
    }

    // Applies above Nyquist conditions to a given MF
    inline void apply_nyquist_conditions(amrex::cMultiFab &field);

};

#include "InflationConfig.impl.hpp"

#endif /* INFLATIONCONFIG */
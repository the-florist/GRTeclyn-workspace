/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

// AMReX includes
#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Reduce.H>

// Doctest includes
#include "doctest.h"
#include "doctestCLIArgs.hpp"

#include <cmath>
#include <vector>

// Test header
#include "InflatonFieldStatisticsTest.hpp"

// Class under test (pulls in InflatonUtils, Potential, StateVariables, ...)
#include "ScalarTensorInit.hpp"

void run_inflaton_field_statistics_test()
{
    int amrex_argc    = doctest::cli_args.argc();
    char **amrex_argv = doctest::cli_args.argv();

    // NOLINTNEXTLINE(bugprone-casting-through-void) // Open MPI triggers this
    amrex::Initialize(amrex_argc, amrex_argv);
    {
        constexpr int num_cells = 32;

        // Populate the parameter table that InflatonParameters and Potential
        // read. Only the tensor sector is exercised: the scalar sector is not
        // yet wired into the state, so it would produce zero perturbations.
        {
            amrex::ParmParse pp;
            pp.addarr("amr.n_cell",
                      std::vector<int>{num_cells, num_cells, num_cells});
            pp.addarr("geometry.prob_extent", std::vector<amrex::Real>{1.0});
            pp.add("scalar_field.G_Newton", 1.0);
        }
        {
            amrex::ParmParse init_pp("init");
            init_pp.add("scalar_init", 0);
            init_pp.add("tensor_init", 1);
            init_pp.add("use_rand", 1);
            init_pp.add("use_window", 0);
            init_pp.add("background_phi", 1.0);
        }
        {
            amrex::ParmParse scalar_field_pp("scalar_field");
            scalar_field_pp.add("scalar_mass", 0.1);
        }

        // Build a single-box 32^3 grid.
        amrex::Box domain{amrex::IntVect::TheZeroVector(),
                          amrex::IntVect{num_cells - 1}};
        amrex::BoxArray box_array{domain};
        amrex::DistributionMapping distribution_mapping{box_array};
        amrex::MFInfo mf_info;
        mf_info.SetArena(amrex::The_Managed_Arena());

        // Zero-initialised state: init() adds perturbations on top, so the
        // state ends up holding the perturbations alone (no background).
        amrex::MultiFab state{box_array, distribution_mapping, NUM_VARS, 0,
                              mf_info};
        state.setVal(0.0);

        // Generate the stochastic initial perturbations.
        ScalarTensorInit initialiser;
        initialiser.init(state);

        // Accumulate the first four raw moments of a tensor perturbation
        // component (h_11) over the whole grid in a single reduction.
        amrex::ReduceOps<amrex::ReduceOpSum, amrex::ReduceOpSum,
                         amrex::ReduceOpSum, amrex::ReduceOpSum>
            reduce_ops;
        amrex::ReduceData<amrex::Real, amrex::Real, amrex::Real, amrex::Real>
            reduce_data(reduce_ops);
        using ReduceTuple = typename decltype(reduce_data)::Type;

        const auto &state_arrays = state.const_arrays();
        reduce_ops.eval(
            state, amrex::IntVect(0), reduce_data,
            [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k) -> ReduceTuple
            {
                const amrex::Real x = state_arrays[box_no](i, j, k, c_h11);
                return {x, x * x, x * x * x, x * x * x * x};
            });

        // reduce_data.value() is local to this rank. With a single-box grid
        // under MPI, only one rank owns the box, so every other rank's local
        // sums are exactly zero (the Sum identity); reduce across ranks to
        // get the true totals everywhere.
        const ReduceTuple moments = reduce_data.value();
        amrex::Real sums[4] = {amrex::get<0>(moments), amrex::get<1>(moments),
                               amrex::get<2>(moments), amrex::get<3>(moments)};
        amrex::ParallelDescriptor::ReduceRealSum(sums, 4);
        const amrex::Real sum_1 = sums[0];
        const amrex::Real sum_2 = sums[1];
        const amrex::Real sum_3 = sums[2];
        const amrex::Real sum_4 = sums[3];

        const amrex::Real num_samples = std::pow(num_cells, 3);

        // Convert raw moments to the mean and central moments.
        const amrex::Real mean      = sum_1 / num_samples;
        const amrex::Real central_2 = sum_2 / num_samples - mean * mean;
        const amrex::Real central_3 = sum_3 / num_samples -
                                      3. * mean * sum_2 / num_samples +
                                      2. * mean * mean * mean;
        const amrex::Real central_4 = sum_4 / num_samples -
                                      4. * mean * sum_3 / num_samples +
                                      6. * mean * mean * sum_2 / num_samples -
                                      3. * mean * mean * mean * mean;

        const amrex::Real standard_deviation = std::sqrt(central_2);
        const amrex::Real skewness = central_3 / std::pow(central_2, 1.5);
        const amrex::Real kurtosis = central_4 / (central_2 * central_2);

        INFO("mean = " << mean << ", stddev = " << standard_deviation
                       << ", skewness = " << skewness
                       << ", kurtosis = " << kurtosis);

        // The perturbations must be present (non-trivial variance).
        CHECK(central_2 > 0.);

        // Mean zero: the k=0 (DC) Fourier mode is zeroed, so the configuration
        // -space mean vanishes to Fourier-transform round-off.
        CHECK(std::abs(mean) < 1.e-6 * standard_deviation);

        // Gaussian: vanishing skewness and an excess kurtosis of zero (i.e. a
        // kurtosis of 3). Tolerances allow for finite-sample scatter over the
        // 32^3 correlated samples.
        constexpr amrex::Real skewness_tolerance = 0.2;
        constexpr amrex::Real kurtosis_tolerance = 0.5;
        CHECK(std::abs(skewness) < skewness_tolerance);
        CHECK(std::abs(kurtosis - 3.) < kurtosis_tolerance);
    }
    amrex::Finalize();
}

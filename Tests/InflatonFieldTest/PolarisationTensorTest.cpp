/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

// AMReX includes
#include <AMReX.H>
#include <AMReX_ParmParse.H>

// Doctest includes
#include "doctest.h"
#include "doctestCLIArgs.hpp"

#include <cmath>
#include <vector>

// Test header
#include "PolarisationTensorTest.hpp"

// Class under test
#include "InflatonUtils.hpp"

void run_polarisation_tensor_test()
{
    int amrex_argc    = doctest::cli_args.argc();
    char **amrex_argv = doctest::cli_args.argv();

    // NOLINTNEXTLINE(bugprone-casting-through-void) // Open MPI triggers this
    amrex::Initialize(amrex_argc, amrex_argv);
    {
        constexpr int num_cells = 8;

        // Populate the parameter table that InflatonParameters and Potential
        // read.
        {
            amrex::ParmParse pp;
            pp.addarr("amr.n_cell",
                      std::vector<int>{num_cells, num_cells, num_cells});
            pp.addarr("geometry.prob_extent", std::vector<amrex::Real>{1.0});
        }
        {
            amrex::ParmParse init_pp("init");
            init_pp.add("background_phi", 1.0);
        }
        {
            amrex::ParmParse scalar_field_pp("scalar_field");
            scalar_field_pp.add("scalar_mass", 1.0);
        }

        const InflatonUtils utils;
        const amrex::Real tol = InflatonUtils::tolerance;

        // Sweep the same half-space of Fourier modes used to build the
        // stochastic fields (i in [0, N/2], j, k in [0, N)), skipping the
        // zero mode where the basis vectors are undefined by construction.
        for (int i = 0; i <= num_cells / 2; ++i)
        {
            for (int j = 0; j < num_cells; ++j)
            {
                for (int k = 0; k < num_cells; ++k)
                {
                    const amrex::IntVect iv{i, j, k};
                    if (iv == amrex::IntVect{0, 0, 0})
                    {
                        continue;
                    }

                    INFO("iv = " << iv);

                    // The basis vectors mhat, nhat must be orthonormal.
                    const auto [mhat, nhat] = utils.calculate_basis_vectors(iv);
                    amrex::Real dot_mm = 0., dot_nn = 0., dot_mn = 0.;
                    for (int l = 0; l < 3; ++l)
                    {
                        dot_mm += mhat[l] * mhat[l];
                        dot_nn += nhat[l] * nhat[l];
                        dot_mn += mhat[l] * nhat[l];
                    }
                    CHECK(std::abs(dot_mm - 1.) < tol);
                    CHECK(std::abs(dot_nn - 1.) < tol);
                    CHECK(std::abs(dot_mn) < tol);

                    // The plus/cross polarisation tensors built from mhat,
                    // nhat must be orthonormal (up to the factor of 2 in the
                    // eplus:eplus and ecross:ecross contractions) and,
                    // individually, trace-free.
                    const auto [eplus, ecross] =
                        utils.calculate_polarisation_tensors(iv);
                    amrex::Real dot_pp = 0., dot_cc = 0., dot_pc = 0.;
                    amrex::Real trace_p = 0., trace_c = 0.;
                    for (int l = 0; l < 3; ++l)
                    {
                        trace_p += eplus(l, l);
                        trace_c += ecross(l, l);
                        for (int p = 0; p < 3; ++p)
                        {
                            dot_pp += eplus(l, p) * eplus(l, p);
                            dot_cc += ecross(l, p) * ecross(l, p);
                            dot_pc += eplus(l, p) * ecross(l, p);
                        }
                    }
                    CHECK(std::abs(dot_pp / 2. - 1.) < tol);
                    CHECK(std::abs(dot_cc / 2. - 1.) < tol);
                    CHECK(std::abs(dot_pc) < tol);
                    CHECK(std::abs(trace_p) < tol);
                    CHECK(std::abs(trace_c) < tol);
                }
            }
        }
    }
    amrex::Finalize();
}

/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

// AMReX includes
#include <AMReX.H>
#include <AMReX_Math.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>

// Doctest includes
#include "doctest.h"
#include "doctestCLIArgs.hpp"

#include <array>
#include <cmath>

// Test header
#include "InflatonFieldBackgroundTest.hpp"

// GRTeclyn includes. Potential.hpp must precede InitialBackgroundData.hpp,
// which uses (but does not include) Potential.
#include "CCZ4RHSWithMatter.hpp"
#include "FourthOrderDerivatives.hpp"
#include "IntegratedMovingPunctureGauge.hpp"
#include "Potential.hpp"
#include "ScalarField.hpp"

#include "InitialBackgroundData.hpp" // pulls in StateVariables.hpp

namespace
{
// A homogeneous state is fully described by the values of the NUM_VARS
// components in a single cell.
using StateVector = std::array<amrex::Real, NUM_VARS>;
} // namespace

void run_inflaton_field_background_test()
{
    int amrex_argc    = doctest::cli_args.argc();
    char **amrex_argv = doctest::cli_args.argv();

    // NOLINTNEXTLINE(bugprone-casting-through-void) // Open MPI triggers this
    amrex::Initialize(amrex_argc, amrex_argv);
    {
        const amrex::Real pi = amrex::Math::pi<amrex::Real>();

        // Quadratic inflation in the slow-roll regime: super-Planckian field
        // amplitude and a light mass.
        constexpr amrex::Real g_newton       = 1.0;
        constexpr amrex::Real scalar_mass    = 0.01;
        constexpr amrex::Real background_phi = 5.0;

        // Start on the slow-roll attractor so there is no initial transient:
        // 3 H phidot = -V' gives a constant slow-roll velocity for a quadratic
        // potential.
        const amrex::Real slow_roll_phidot =
            -scalar_mass / std::sqrt(12. * pi * g_newton);

        // Populate the parameter table read by InitialBackgroundData, the
        // Potential, the ScalarField matter and the CCZ4 RHS / gauge.
        {
            amrex::ParmParse pp;
            pp.add("scalar_field.G_Newton", g_newton);
            pp.add("evolution.sigma", 0.0);
        }
        {
            amrex::ParmParse init_pp("init");
            init_pp.add("background_phi", background_phi);
            init_pp.add("background_dphi", slow_roll_phidot);
        }
        {
            amrex::ParmParse scalar_field_pp("scalar_field");
            scalar_field_pp.add("scalar_mass", scalar_mass);
        }
        {
            // BSSN formulation with all damping/dissipation off.
            amrex::ParmParse ccz4_pp("ccz4");
            ccz4_pp.add("formulation", 1);
            ccz4_pp.add("kappa1", 0.0);
            ccz4_pp.add("kappa2", 0.0);
            ccz4_pp.add("kappa3", 0.0);
            ccz4_pp.add("covariantZ4", 0);
            ccz4_pp.add("min_chi", 1.e-30);
            ccz4_pp.add("min_lapse", 1.e-30);
        }
        {
            // Frozen gauge: the lapse stays at one, so coordinate time is
            // proper (cosmic) time and the shift stays zero.
            amrex::ParmParse gauge_pp("gauge");
            gauge_pp.add("lapse_advec_coeff", 0.0);
            gauge_pp.add("lapse_coeff", 0.0);
            gauge_pp.add("lapse_power", 1.0);
            gauge_pp.add("shift_advec_coeff", 0.0);
            gauge_pp.add("shift_Gamma_coeff", 0.0);
            gauge_pp.add("eta", 0.0);
        }

        // A small cubic grid. Every cell holds the same homogeneous data, so
        // filling the ghost cells to the same value makes every spatial
        // derivative vanish and the CCZ4 + matter RHS reduces to the FLRW ODEs.
        constexpr int num_cells  = 4;
        constexpr int num_ghosts = 3;
        constexpr amrex::Real dx = 1.0; // irrelevant: derivatives vanish

        amrex::Box domain{amrex::IntVect::TheZeroVector(),
                          amrex::IntVect{num_cells - 1}};
        amrex::BoxArray box_array{domain};
        amrex::DistributionMapping distribution_mapping{box_array};
        amrex::MFInfo mf_info;
        mf_info.SetArena(amrex::The_Managed_Arena());

        amrex::MultiFab soln_mf{box_array, distribution_mapping, NUM_VARS,
                                num_ghosts, mf_info};
        amrex::MultiFab rhs_mf{box_array, distribution_mapping, NUM_VARS,
                               num_ghosts, mf_info};

        // Set the flat, unit-lapse FLRW background via InitialBackgroundData
        // and read it back into a single state vector.
        StateVector state{};
        {
            soln_mf.setVal(0.0);
            InitialBackgroundData background;
            const auto &soln_arrays = soln_mf.arrays();
            amrex::ParallelFor(
                soln_mf, soln_mf.nGrowVect(),
                [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k) noexcept
                { background.compute(i, j, k, soln_arrays[box_no]); });
            amrex::Gpu::streamSynchronize();

            for (amrex::MFIter mfi(soln_mf); mfi.isValid(); ++mfi)
            {
                const auto array = soln_mf.const_array(mfi);
                for (int n = 0; n < NUM_VARS; ++n)
                {
                    state[n] = array(0, 0, 0, n);
                }
                break;
            }
        }

        // RHS operator: fill the grid uniformly with the given homogeneous
        // state and evaluate the production CCZ4 + matter RHS at one cell.
        const CCZ4RHSWithMatter<ScalarField<Potential, FourthOrderDerivatives>,
                                FourthOrderDerivatives>
            ccz4_rhs{dx};
        const IntegratedMovingPunctureGauge<FourthOrderDerivatives> gauge{dx};

        auto evaluate_rhs = [&](const StateVector &value) -> StateVector
        {
            for (int n = 0; n < NUM_VARS; ++n)
            {
                soln_mf.setVal(value[n], n, 1, num_ghosts);
            }
            rhs_mf.setVal(0.0);

            const auto &soln_arrays = soln_mf.const_arrays();
            const auto &rhs_arrays  = rhs_mf.arrays();
            amrex::ParallelFor(
                rhs_mf,
                [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k)
                {
                    ccz4_rhs.compute_chi_and_h_ij(i, j, k, rhs_arrays[box_no],
                                                  soln_arrays[box_no]);
                    ccz4_rhs.compute_A_ij_and_Theta_and_Gamma(
                        i, j, k, rhs_arrays[box_no], soln_arrays[box_no]);
                    ccz4_rhs.add_emtensor_rhs(i, j, k, rhs_arrays[box_no],
                                              soln_arrays[box_no]);
                    gauge.calculate_rhs(i, j, k, rhs_arrays[box_no],
                                        soln_arrays[box_no]);
                    ccz4_rhs.add_matter_rhs(i, j, k, rhs_arrays[box_no],
                                            soln_arrays[box_no]);
                });
            amrex::Gpu::streamSynchronize();

            StateVector rhs{};
            for (amrex::MFIter mfi(rhs_mf); mfi.isValid(); ++mfi)
            {
                const auto array = rhs_mf.const_array(mfi);
                for (int n = 0; n < NUM_VARS; ++n)
                {
                    rhs[n] = array(0, 0, 0, n);
                }
                break;
            }
            return rhs;
        };

        // Explicit RK4 integration for a few steps.
        constexpr int num_steps  = 40;
        constexpr amrex::Real dt = 0.25;

        auto axpy = [](const StateVector &v, const StateVector &k,
                       const amrex::Real factor) -> StateVector
        {
            StateVector out{};
            for (int n = 0; n < NUM_VARS; ++n)
            {
                out[n] = v[n] + factor * k[n];
            }
            return out;
        };

        for (int step = 0; step < num_steps; ++step)
        {
            const StateVector k1 = evaluate_rhs(state);
            const StateVector k2 = evaluate_rhs(axpy(state, k1, 0.5 * dt));
            const StateVector k3 = evaluate_rhs(axpy(state, k2, 0.5 * dt));
            const StateVector k4 = evaluate_rhs(axpy(state, k3, dt));
            for (int n = 0; n < NUM_VARS; ++n)
            {
                state[n] += dt / 6. * (k1[n] + 2. * k2[n] + 2. * k3[n] + k4[n]);
            }
        }

        // Evolved values from the code.
        const amrex::Real phi_code    = state[c_phi];
        const amrex::Real pi_code     = state[c_Pi];
        const amrex::Real k_code      = state[c_K];
        const amrex::Real hubble_code = -k_code / 3.;
        const amrex::Real potential_code =
            0.5 * scalar_mass * scalar_mass * phi_code * phi_code;

        // Analytic slow-roll solution at the final time.
        const amrex::Real final_time = num_steps * dt;
        const amrex::Real phi_slow_roll =
            background_phi + slow_roll_phidot * final_time;
        const amrex::Real potential_slow_roll =
            0.5 * scalar_mass * scalar_mass * phi_slow_roll * phi_slow_roll;
        const amrex::Real hubble_slow_roll =
            std::sqrt(8. * pi * g_newton / 3. * potential_slow_roll);

        INFO("phi: code = " << phi_code << ", slow-roll = " << phi_slow_roll
                            << "; H: code = " << hubble_code
                            << ", slow-roll = " << hubble_slow_roll);

        // The background must actually have evolved.
        CHECK(std::abs(phi_code - background_phi) >
              0.5 * std::abs(slow_roll_phidot * final_time));

        // The exact Friedmann (Hamiltonian) constraint,
        // H^2 = (8 pi G / 3)(1/2 Pi^2 + V), must hold throughout the evolution.
        const amrex::Real friedmann_rhs =
            8. * pi * g_newton / 3. *
            (0.5 * pi_code * pi_code + potential_code);
        CHECK(std::abs(hubble_code * hubble_code - friedmann_rhs) <
              1.e-3 * hubble_code * hubble_code);

        // Slow-roll: the kinetic energy is subdominant to the potential.
        CHECK(0.5 * pi_code * pi_code < 1.e-2 * potential_code);

        // The solution follows the analytic slow-roll trajectory.
        CHECK(std::abs(phi_code - phi_slow_roll) <
              2.e-2 * std::abs(phi_slow_roll));
        CHECK(std::abs(hubble_code - hubble_slow_roll) <
              2.e-2 * std::abs(hubble_slow_roll));
    }
    amrex::Finalize();
}

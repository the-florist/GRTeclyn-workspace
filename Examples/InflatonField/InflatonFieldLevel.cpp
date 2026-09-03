/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#include "InflatonFieldLevel.hpp"

#include "AlgebraicConstraintsEnforcer.hpp"
#include "CCZ4RHSWithMatter.hpp"
#include "ConstraintsWithMatter.hpp"
#include "EMTensor.hpp"
#include "FixedGridsTagger.hpp"
#include "FourthOrderDerivatives.hpp"
#include "GammaCalculator.hpp"
#include "IntegratedMovingPunctureGauge.hpp"
#include "PositiveChiAndLapse.hpp"
#include "SixthOrderDerivatives.hpp"
#include "StateTypes.hpp"

#include <algorithm>
#include <type_traits>

// // Problem specific includes
#include "DerivedVariables.hpp"
#include "InitialBackgroundData.hpp"
#include "Potential.hpp"
#include "ScalarField.hpp"
#include "ScalarTensorInit.hpp"

using InflatonFieldConstraints =
    ConstraintsWithMatter<InflatonFieldLevel::ScalarFieldWithPotential<>>;

void InflatonFieldLevel::variableSetUp()
{
    BL_PROFILE("InflatonFieldLevel::variableSetUp()");
    state_variable_set_up();
    InflatonFieldConstraints::set_up(state_index);
    DerivedVariables::set_up(state_index);

    // Ensure that if the user requests any of the derived
    // variables in DerivedVariables, they are NOT using AMR
    // at the same time.
    GRParmParse pp("amr");
    amrex::Vector<std::string> derive_plot_vars;
    int max_level = -1;
    pp.queryarr("derive_plot_vars", derive_plot_vars);
    pp.get("max_level", max_level);

    const bool all_requested =
        std::find(derive_plot_vars.begin(), derive_plot_vars.end(), "ALL") !=
        derive_plot_vars.end();

    bool derive_requested = false;
    for (const auto &name : DerivedVariables::var_names)
    {
        if (std::find(derive_plot_vars.begin(), derive_plot_vars.end(), name) !=
            derive_plot_vars.end())
        {
            derive_requested = true;
            break;
        }
    }

    if ((derive_requested || all_requested) && max_level != 0)
    {
        amrex::Error("InflatonFieldLevel::variableSetUp, "
                     "DerivedVariables cannot be used with AMR");
    }
}

// Things to do at each advance step, after the RK4 is calculated
void InflatonFieldLevel::specific_advance()
{
    BL_PROFILE("InflatonFieldLevel::specific_advance()");

    amrex::MultiFab &state_new = get_new_data(state_index);
    const auto &state_arrays   = state_new.arrays();

    const AlgebraicConstraintsEnforcer algebraic_constraints_enforcer;
    const PositiveChiAndLapse positive_chi_and_lapse;

    amrex::ParallelFor(
        state_new,
        [=] AMREX_GPU_DEVICE(int box_no, int ix, int iy, int iz)
        {
            algebraic_constraints_enforcer(ix, iy, iz, state_arrays[box_no]);
            positive_chi_and_lapse(ix, iy, iz, state_arrays[box_no]);
        });
    amrex::Gpu::streamSynchronize();
}

// Initial data for field and metric variables
void InflatonFieldLevel::initData()
{
    BL_PROFILE("InflatonFieldLevel::initData()");
    if (get_gr_amr_ptr()->Verbose() > 0)
    {
        amrex::Print() << "InflatonFieldLevel::initData " << Level() << "\n";
    }

    amrex::MultiFab &state_new = get_new_data(state_index);
    const auto &state_arrays   = state_new.arrays();

    InitialBackgroundData FLRW_background;
    amrex::ParallelFor(
        state_new, state_new.nGrowVect(),
        [=] AMREX_GPU_DEVICE(int box_ind, int i, int j, int k) noexcept
        {
            amrex::CellData<amrex::Real> cell =
                state_arrays[box_ind].cellData(i, j, k);
            for (int n = 0; n < cell.nComp(); ++n)
            {
                cell[n] = 0.;
            }

            FLRW_background.compute(i, j, k, state_arrays[box_ind]);
        });

    amrex::Gpu::streamSynchronize();

    ScalarTensorInit random_field_initialiser;
    random_field_initialiser.init(state_new);

    if (m_evolution_spatial_derivative_order == 4)
    {
        const GammaCalculator<FourthOrderDerivatives> gamma_calculator(
            Geom().CellSize(0));
        const IntegratedMovingPunctureGauge<FourthOrderDerivatives>
            integrated_moving_puncture_gauge(Geom().CellSize(0));
        amrex::ParallelFor(
            state_new,
            [=] AMREX_GPU_DEVICE(int box_no, int ix, int iy, int iz)
            {
                gamma_calculator(ix, iy, iz, state_arrays[box_no]);
                integrated_moving_puncture_gauge.set_initial_B_to_Gamma(
                    ix, iy, iz, state_arrays[box_no]);
            });
    }
    else if (m_evolution_spatial_derivative_order == 6)
    {
        const GammaCalculator<SixthOrderDerivatives> gamma_calculator(
            Geom().CellSize(0));
        const IntegratedMovingPunctureGauge<SixthOrderDerivatives>
            integrated_moving_puncture_gauge(Geom().CellSize(0));
        amrex::ParallelFor(
            state_new,
            [=] AMREX_GPU_DEVICE(int box_no, int ix, int iy, int iz)
            {
                gamma_calculator(ix, iy, iz, state_arrays[box_no]);
                integrated_moving_puncture_gauge.set_initial_B_to_Gamma(
                    ix, iy, iz, state_arrays[box_no]);
            });
    }

    amrex::Gpu::streamSynchronize();
}

// NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
void InflatonFieldLevel::specific_eval_rhs(amrex::MultiFab &a_soln,
                                           amrex::MultiFab &a_rhs,
                                           const amrex::Real /*a_time*/)
{
    BL_PROFILE("InflatonFieldLevel::specific_eval_rhs()");

    const auto &soln_arrays       = a_soln.arrays();
    const auto &const_soln_arrays = a_soln.const_arrays();
    const auto &rhs_arrays        = a_rhs.arrays();

    const AlgebraicConstraintsEnforcer algebraic_constraints_enforcer;
    const PositiveChiAndLapse positive_chi_and_lapse;

    amrex::ParallelFor(
        a_soln, a_soln.nGrowVect(),
        [=] AMREX_GPU_DEVICE(int box_no, int ix, int iy, int iz)
        {
            algebraic_constraints_enforcer(ix, iy, iz, soln_arrays[box_no]);
            positive_chi_and_lapse(ix, iy, iz, soln_arrays[box_no]);
        });

    if (m_evolution_spatial_derivative_order != 4 &&
        m_evolution_spatial_derivative_order != 6)
    {
        amrex::Abort("spatial_derivative_order must be 4 or 6");
    }

    // NOLINTBEGIN(bugprone-branch-clone)
    if (m_evolution_spatial_derivative_order == 4)
    {
        const CCZ4RHSWithMatter<
            ScalarFieldWithPotential<FourthOrderDerivatives>,
            FourthOrderDerivatives>
            ccz4_rhs(Geom().CellSize(0));
        const IntegratedMovingPunctureGauge<FourthOrderDerivatives>
            integrated_moving_puncture_gauge(Geom().CellSize(0));

        amrex::ParallelFor(
            a_rhs,
            [=] AMREX_GPU_DEVICE(int box_no, int ix, int iy, int iz)
            {
                ccz4_rhs.compute_chi_and_h_ij(ix, iy, iz, rhs_arrays[box_no],
                                              const_soln_arrays[box_no]);
            });

        // NOLINTBEGIN(bugprone-easily-swappable-parameters)
        amrex::ParallelFor(
            a_rhs,
            [=] AMREX_GPU_DEVICE(int box_no, int ix, int iy, int iz)
            {
                ccz4_rhs.compute_A_ij_and_Theta_and_Gamma(
                    ix, iy, iz, rhs_arrays[box_no], const_soln_arrays[box_no]);
            });
        // NOLINTEND(bugprone-easily-swappable-parameters)

        amrex::ParallelFor(
            a_rhs,
            [=] AMREX_GPU_DEVICE(int box_no, int ix, int iy, int iz)
            {
                ccz4_rhs.add_emtensor_rhs(ix, iy, iz, rhs_arrays[box_no],
                                          const_soln_arrays[box_no]);
                integrated_moving_puncture_gauge.calculate_rhs(
                    ix, iy, iz, rhs_arrays[box_no], const_soln_arrays[box_no]);
                ccz4_rhs.add_matter_rhs(ix, iy, iz, rhs_arrays[box_no],
                                        const_soln_arrays[box_no]);
                ccz4_rhs.apply_dissipation(ix, iy, iz, rhs_arrays[box_no],
                                           const_soln_arrays[box_no]);
            });
    }
    else if (m_evolution_spatial_derivative_order == 6)
    {
        const CCZ4RHSWithMatter<ScalarFieldWithPotential<SixthOrderDerivatives>,
                                SixthOrderDerivatives>
            ccz4_rhs(Geom().CellSize(0));
        const IntegratedMovingPunctureGauge<SixthOrderDerivatives>
            integrated_moving_puncture_gauge(Geom().CellSize(0));

        amrex::ParallelFor(
            a_rhs,
            [=] AMREX_GPU_DEVICE(int box_no, int ix, int iy, int iz)
            {
                ccz4_rhs.compute_chi_and_h_ij(ix, iy, iz, rhs_arrays[box_no],
                                              const_soln_arrays[box_no]);
            });

        // NOLINTBEGIN(bugprone-easily-swappable-parameters)
        amrex::ParallelFor(
            a_rhs,
            [=] AMREX_GPU_DEVICE(int box_no, int ix, int iy, int iz)
            {
                ccz4_rhs.compute_A_ij_and_Theta_and_Gamma(
                    ix, iy, iz, rhs_arrays[box_no], const_soln_arrays[box_no]);
            });
        // NOLINTEND(bugprone-easily-swappable-parameters)

        amrex::ParallelFor(
            a_rhs,
            [=] AMREX_GPU_DEVICE(int box_no, int ix, int iy, int iz)
            {
                ccz4_rhs.add_emtensor_rhs(ix, iy, iz, rhs_arrays[box_no],
                                          const_soln_arrays[box_no]);
                integrated_moving_puncture_gauge.calculate_rhs(
                    ix, iy, iz, rhs_arrays[box_no], const_soln_arrays[box_no]);
                ccz4_rhs.add_matter_rhs(ix, iy, iz, rhs_arrays[box_no],
                                        const_soln_arrays[box_no]);
                ccz4_rhs.apply_dissipation(ix, iy, iz, rhs_arrays[box_no],
                                           const_soln_arrays[box_no]);
            });
    }
    // NOLINTEND(bugprone-branch-clone)

    amrex::Gpu::streamSynchronize();
}

void InflatonFieldLevel::specific_update_ode(amrex::MultiFab &a_soln)
{
    BL_PROFILE("InflatonFieldLevel::specific_update_ode()");

    const auto &soln_arrays = a_soln.arrays();
    const AlgebraicConstraintsEnforcer algebraic_constraints_enforcer;

    amrex::ParallelFor(
        a_soln, amrex::IntVect(0),
        [=] AMREX_GPU_DEVICE(int box_no, int ix, int iy, int iz)
        { algebraic_constraints_enforcer(ix, iy, iz, soln_arrays[box_no]); });
    amrex::Gpu::streamSynchronize();
}

void InflatonFieldLevel::tag_cells(amrex::TagBoxArray &a_tag_box_array,
                                   const amrex::Real /*a_regrid_threshold*/)
{
    BL_PROFILE("InflatonFieldLevel::tag_cells()");

    const auto &tag_arrays = a_tag_box_array.arrays();

    const FixedGridsTagger tagger(Geom().CellSize(0), Level());

    amrex::ParallelFor(a_tag_box_array,
                       [=] AMREX_GPU_DEVICE(int box_no, int ix, int iy, int iz)
                       { tagger(ix, iy, iz, tag_arrays[box_no]); });
    amrex::Gpu::streamSynchronize();
}

void InflatonFieldLevel::specific_post_timestep()
{
    BL_PROFILE("InflatonFieldLevel::specific_post_timestep()");
}

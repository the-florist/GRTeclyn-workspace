/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#include "ScalarFieldLevel.hpp"

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

#include <type_traits>

// // Problem specific includes
#include "InitialBackgroundData.hpp"
#include "Potential.hpp"
#include "ScalarField.hpp"
#include "RandomFieldInit.hpp"
#include "InflationExtraction.hpp"

using ScalarFieldConstraints =
    ConstraintsWithMatter<ScalarFieldLevel::ScalarFieldWithPotential<>>;

void ScalarFieldLevel::variableSetUp()
{
    BL_PROFILE("ScalarFieldLevel::variableSetUp()");
    state_variable_set_up();
    ScalarFieldConstraints::set_up(state_index);
    InflationExtraction::set_up(state_index);
}

// Things to do at each advance step, after the RK4 is calculated
void ScalarFieldLevel::specific_advance()
{
    BL_PROFILE("ScalarFieldLevel::specific_advance()");

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
void ScalarFieldLevel::initData()
{
    BL_PROFILE("ScalarFieldLevel::initData()");
    if (get_gr_amr_ptr()->Verbose() > 0)
    {
        amrex::Print() << "ScalarFieldLevel::initData " << Level() << "\n";
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

    RandomFieldInit random_field_initialiser;
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
void ScalarFieldLevel::specific_eval_rhs(amrex::MultiFab &a_soln,
                                         amrex::MultiFab &a_rhs,
                                         const amrex::Real /*a_time*/)
{
    BL_PROFILE("ScalarFieldLevel::specific_eval_rhs()");

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

void ScalarFieldLevel::specific_update_ode(amrex::MultiFab &a_soln)
{
    BL_PROFILE("ScalarFieldLevel::specific_update_ode()");

    const auto &soln_arrays = a_soln.arrays();
    const AlgebraicConstraintsEnforcer algebraic_constraints_enforcer;

    amrex::ParallelFor(
        a_soln, amrex::IntVect(0),
        [=] AMREX_GPU_DEVICE(int box_no, int ix, int iy, int iz)
        { algebraic_constraints_enforcer(ix, iy, iz, soln_arrays[box_no]); });
    amrex::Gpu::streamSynchronize();
}

void ScalarFieldLevel::tag_cells(amrex::TagBoxArray &a_tag_box_array,
                                 const amrex::Real /*a_regrid_threshold*/)
{
    BL_PROFILE("ScalarFieldLevel::tag_cells()");

    const auto &tag_arrays = a_tag_box_array.arrays();

    const FixedGridsTagger tagger(Geom().CellSize(0), Level());

    amrex::ParallelFor(a_tag_box_array,
                       [=] AMREX_GPU_DEVICE(int box_no, int ix, int iy, int iz)
                       { tagger(ix, iy, iz, tag_arrays[box_no]); });
    amrex::Gpu::streamSynchronize();
}

void ScalarFieldLevel::specific_post_timestep()
{
	BL_PROFILE("ScalarFieldLevel::specific_post_timestep()");
    if (Level() == 0)
    {
        GRParmParse pp;
        const amrex::Real time         = get_state_data(state_index).curTime();
        const amrex::Real dt           = get_gr_amr_ptr()->dtLevel(0);
        const amrex::Real restart_time = get_gr_amr_ptr()->get_restart_time();
        const bool first_step          = (time <= dt);

        amrex::MultiFab &state_new = get_new_data(state_index);

        std::string output_path, data_path, print_path;
        pp.get("output_path", data_path);
        pp.get("data_path", data_path);
        print_path = output_path + data_path;
        
        int N;
        pp.get("N", N);

        const int vol = std::pow(N, 3.); // (!!) unitless volume
        const double phi_avg = state_new.sum(c_phi)/vol;
        const double Pi_avg = state_new.sum(c_Pi)/vol;
        const double chi_avg = state_new.sum(c_chi)/vol;
        const double scale_fact_avg = 1./sqrt(chi_avg);
        const double Hubble_fact_avg = -state_new.sum(c_K)/vol/3.;
        const double lapse_avg = state_new.sum(c_lapse)/vol;

        const amrex::BoxArray& ba = state_new.boxArray();
        const amrex::DistributionMapping& dm = state_new.DistributionMap();
        int ngrow = state_new.nGrow();

        SmallDataIO means_file(print_path+"means-file", dt, time, restart_time, SmallDataIO::APPEND, first_step, ".dat");
        means_file.remove_duplicate_time_data(); // removes any duplicate data from previous run (for checkpointing)

        if(first_step) 
        {
            means_file.write_header_line({"PhiMean","PiMean","ScaleFactMean","HubbleMean","LapseMean"});
        }
        const std::vector<amrex::Real> means_data = 
            {phi_avg, Pi_avg, scale_fact_avg, Hubble_fact_avg, lapse_avg};
        means_file.write_time_data_line(means_data);

        // Extract the spectra and field statistics
        InflationExtraction inflation_extractor_engine;
        inflation_extractor_engine.set_print_params(print_path, time, dt, 
                                                    restart_time);
        inflation_extractor_engine.extract(state_new);

        // Make a file object for constraint statistics
        SmallDataIO constrs_file(print_path+"constraint-statistics", dt, 
                                time, restart_time, SmallDataIO::APPEND, first_step, ".dat");
        constrs_file.remove_duplicate_time_data();
        
        // Find the constraints and put them in a MF
        int num = 0;
        const std::list<amrex::DeriveRec>& dlist = derive_lst.dlist();
        for (auto const& var: dlist)
        {
            if(var.name() == "constraints") { num = var.numDerive(); }
        }
        if(first_step) { std::cout << "Num derive vars: " << num << "\n"; }

        amrex::MultiFab constr_alias(ba, dm, num, ngrow, amrex::MFInfo(), Factory());
        constr_alias.setVal(0.0);
        derive("constraints", time, constr_alias, 0);

        // Print statistics on the abs constraint terms
        amrex::Vector<int> moments{1,2};
        inflation_extractor_engine.print_moment(constr_alias, Constraints::var_names, 
                                            moments, constrs_file, first_step);
        }
    }

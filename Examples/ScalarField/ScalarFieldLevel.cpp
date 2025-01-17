/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

// General includes common to most GR problems
#include "ScalarFieldLevel.hpp"
#include "PositiveChiAndAlpha.hpp"
#include "TraceARemoval.hpp"
// //#include "SixthOrderDerivatives.hpp"

// // For RHS update
#include "MatterCCZ4RHS.hpp"

// // For constraints calculation
#include "Constraints.hpp"
#include "MatterWeyl4.hpp"
#include "NewMatterConstraints.hpp"

// For tagging cells
#include "ChiExtractionTaggingCriterion.hpp"
#include "FixedGridsTaggingCriterion.hpp"

// // Problem specific includes
#include "InitialBackgroundData.hpp"
#include "InitialScalarData.hpp"
#include "Potential.hpp"
#include "ScalarField.hpp"
#include "RandomField.hpp"

void ScalarFieldLevel::variableSetUp()
{
    BL_PROFILE("ScalarFieldLevel::variableSetUp()");

    // Set up the state variables
    stateVariableSetUp();

    const int nghost = simParams().num_ghosts;

    // // Add the constraints to the derive list
    derive_lst.add(
        "constraints", amrex::IndexType::TheCellType(),
        static_cast<int>(Constraints::var_names.size()), Constraints::var_names,
        amrex::DeriveFuncFab(), // null function because we won't use
                                // it.
        [=](const amrex::Box &box) { return amrex::grow(box, nghost); },
        &amrex::cell_quartic_interp);

    // We only need the non-gauge CCZ4 variables to calculate the constraints
    derive_lst.addComponent("constraints", desc_lst, State_Type, 0, c_lapse);

    // Add Weyl4 to the derive list
    derive_lst.add(
        "Weyl4", amrex::IndexType::TheCellType(),
        static_cast<int>(Weyl4::var_names.size()), Weyl4::var_names,
        amrex::DeriveFuncFab(), // null function because we won't use it
        [=](const amrex::Box &box) { return amrex::grow(box, nghost); },
        &amrex::cell_quartic_interp);

    // We need all of the CCZ4 variables to calculate Weyl4 (except B)
    derive_lst.addComponent("Weyl4", desc_lst, State_Type, 0, c_B1);
}

// Things to do at each advance step, after the RK4 is calculated
void ScalarFieldLevel::specificAdvance()
{
    BL_PROFILE("ScalarFieldLevel::specificAdvance");
    // Enforce trace free A_ij and positive chi and alpha
    amrex::MultiFab &S_new = get_new_data(State_Type);
    const auto &arrs       = S_new.arrays();

    // Enforce the trace free A_ij condition and positive chi and alpha
    amrex::ParallelFor(S_new,
                       [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k)
                       {
                           amrex::CellData<amrex::Real> cell =
                               arrs[box_no].cellData(i, j, k);
                           TraceARemoval()(cell);
                           PositiveChiAndAlpha()(cell);
                       });

    // Check for nan's
    if (simParams().nan_check)
    {
        if (S_new.contains_nan(0, S_new.nComp(), amrex::IntVect(0), true))
        {
            amrex::Abort("NaN in specificAdvance");
        }
    }
}

// Initial data for field and metric variables
void ScalarFieldLevel::initData()
{
    BL_PROFILE("ScalarFieldLevel::initData");
    if (m_verbosity)
        amrex::Print() << "ScalarFieldLevel::initialData " << Level()
                       << std::endl;

    const auto dx = geom.CellSizeArray();
    InitialBackgroundData FLRW_background(simParams().background_params);

    amrex::MultiFab &state  = get_new_data(State_Type);
    auto const &state_array = state.arrays();

    amrex::ParallelFor(
        state, state.nGrowVect(),
        [=] AMREX_GPU_DEVICE(int box_ind, int i, int j, int k) noexcept
        {
            amrex::CellData<amrex::Real> cell =
                state_array[box_ind].cellData(i, j, k);
            for (int n = 0; n < cell.nComp(); ++n)
            {
                cell[n] = 0.;
            }

            FLRW_background.compute(i, j, k, state_array[box_ind]);
        });

    RandomField random_field(simParams().random_field_params, simParams().background_params, "field");
    random_field.init_random_field();

    if (simParams().nan_check)
    {
        if (state.contains_nan(0, state.nComp(), amrex::IntVect(0), true))
        {
            amrex::Abort("NaN in initData");
        }
    }
}

// Things to do in RHS update, at each RK4 step
void ScalarFieldLevel::specificEvalRHS(amrex::MultiFab &a_soln,
                                       amrex::MultiFab &a_rhs,
                                       const double a_time)
{
    BL_PROFILE("ScalarFieldLevel::specificEvalRHS()");

    // amrex::WriteSingleLevelPlotfileHDF5(plot_file_root, a_rhs,
    // 				      var_names, Geom(), 0.0, 0);

    const auto &soln_arrs   = a_soln.arrays();
    const auto &soln_c_arrs = a_soln.const_arrays();
    const auto &rhs_arrs    = a_rhs.arrays();

    // Enforce positive chi and alpha and trace free A
    amrex::ParallelFor(a_soln, a_soln.nGrowVect(),
                       [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k)
                       {
                           amrex::CellData<amrex::Real> cell =
                               soln_arrs[box_no].cellData(i, j, k);
                           TraceARemoval()(cell);
                           PositiveChiAndAlpha()(cell);
                       });

    // Calculate MatterCCZ4 right hand side with matter_t = ScalarField
    Potential potential(simParams().potential_params);
    ScalarFieldWithPotential scalar_field(potential);

    // Calculate CCZ4 right hand side
    if (simParams().max_spatial_derivative_order == 4)
    {
        MatterCCZ4RHS<ScalarFieldWithPotential, MovingPunctureGauge,
                      FourthOrderDerivatives>
            matter_ccz4_rhs(scalar_field, simParams().ccz4_params,
                            Geom().CellSize(0), simParams().sigma,
                            simParams().formulation, simParams().G_Newton);
        amrex::ParallelFor(
            a_rhs,
            [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k) {
                matter_ccz4_rhs.compute(i, j, k, rhs_arrs[box_no],
                                        soln_c_arrs[box_no]);
            });
    }
    else if (simParams().max_spatial_derivative_order == 6)
    {
        amrex::Abort("xxxxx max_spatial_derivative_order == 6 todo");
#if 0
        MatterCCZ4RHS<ScalarFieldWithPotential, MovingPunctureGauge, SixthOrderDerivatives>
	  matter_ccz4_rhs(scalar_field, simParams().ccz4_params, Geom().CellSize(0), simParams().sigma,
			  simParams().formulation, simParams().G_Newton);
        amrex::ParallelFor(a_rhs,
        [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k)
        {
            amrex::CellData<amrex::Real const> state = soln_c_arrs[box_no].cellData(i,j,k);
            amrex::CellData<amrex::Real> rhs = rhs_arrs[box_no].cellData(i,j,k);
            matter_ccz4_rhs.compute(i,j,k,rhs_arrs[box_no], soln_c_arrs[box_no]);
        });
#endif
    }

    if (simParams().nan_check)
    {
        if (a_soln.contains_nan(0, a_soln.nComp(), amrex::IntVect(0), true))
        {

            for (int i = 0; i < NUM_VARS; i++)
            {
                if (a_soln.contains_nan(i, 1, amrex::IntVect(0), true))
                    amrex::Print() << "Nan found in "
                                   << StateVariables::names[i] << std::endl;
            }

            amrex::Abort("NaN in specificUpdateRHS");
        }
    }
}

// Things to do at ODE update, after soln + rhs
void ScalarFieldLevel::specificUpdateODE(amrex::MultiFab &a_soln)
{
    BL_PROFILE("ScalarFieldLevel::specificUpdateODE()");
    // Enforce the trace free A_ij condition
    const auto &soln_arrs = a_soln.arrays();
    amrex::ParallelFor(a_soln, amrex::IntVect(0), // zero ghost cells
                       [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k)
                       {
                           amrex::CellData<amrex::Real> cell =
                               soln_arrs[box_no].cellData(i, j, k);
                           TraceARemoval()(cell);
                       });
}

void ScalarFieldLevel::preTagCells()
{
    // we don't need any ghosts filled for the fixed grids tagging criterion
    // used here so don't fill any
}

void ScalarFieldLevel::errorEst(amrex::TagBoxArray &tagging_criterion,
                                int /*clearval*/, int /*tagval*/,
                                amrex::Real /*time*/, int /*n_error_buf*/,
                                int /*ngrow*/)

{
    BL_PROFILE("ScalarFieldLevel::errorEst()");

    amrex::MultiFab &state_new = get_new_data(State_Type);
    const auto curr_time       = get_state_data(State_Type).curTime();

    const int nghost =
        state_new.nGrow(); // Need ghost cells to compute gradient
    const int ncomp = state_new.nComp();

    // I filled all the ghost cells in case but could also just fill the ones
    // used for tagging We only use chi in the tagging criterion so only fill
    // the ghosts for chi
    FillPatch(*this, state_new, nghost, curr_time, State_Type, 0, ncomp);

    const auto &simpar = simParams();

    const auto &tags           = tagging_criterion.arrays();
    const auto &state_new_arrs = state_new.const_arrays();
    const auto tagval          = amrex::TagBox::SET;

    amrex::Real dx     = Geom().CellSize(0);
    int curr_level     = Level();
    const auto probhi  = Geom().ProbHiArray();
    const auto problo  = Geom().ProbLoArray();
    amrex::Real length = probhi[0] - problo[0];

    const auto test_dx = Geom().CellSizeArray();

    FixedGridsTaggingCriterion tagger(test_dx[0], curr_level, length,
                                      simParams().initial_params.center);
    // ChiExtractionTaggingCriterion tagger(Geom().CellSize(0), Level(),
    //                                      simpar.extraction_params,
    //                                      simpar.activate_extraction);

    amrex::Real threshold = simpar.regrid_thresholds[Level()];
    amrex::ParallelFor(state_new, amrex::IntVect(0),
                       [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k)
                       {
                           amrex::Real criterion =
                               tagger.compute(i, j, k, state_new_arrs[box_no]);

                           // amrex::Real criterion =
                           //     tagger(i, j, k, state_new_arrs[box_no]);

                           if (criterion >= threshold)
                           {
                               tags[box_no](i, j, k) = tagval;
                           }
                       });

    amrex::Gpu::streamSynchronize();
}

void ScalarFieldLevel::derive(const std::string &name, amrex::Real time,
                              amrex::MultiFab &multifab, int dcomp)
{
    BL_PROFILE("ScalarFieldLevel::derive()");

    BL_ASSERT(dcomp < multifab.nComp());

    const int num_ghosts = multifab.nGrow();

    const amrex::DeriveRec *rec = derive_lst.get(name);
    if (rec != nullptr)
    {
        int state_idx, derive_scomp, derive_ncomp;

        // we only have one state so state_idx will be State_Type = 0
        rec->getRange(0, state_idx, derive_scomp, derive_ncomp);

        // work out how many extra ghost cells we need
        const amrex::BoxArray &src_ba = state[state_idx].boxArray();

        int num_extra_ghosts = num_ghosts;
        {
            amrex::Box box0   = src_ba[0];
            amrex::Box box1   = rec->boxMap()(box0);
            num_extra_ghosts += box0.smallEnd(0) - box1.smallEnd(0);
        }

        // Make a Multifab with enough extra ghosts to calculated derived
        // quantity. For now use NUM_VARS in case the enum mapping loads more
        // vars than is actually needed
        amrex::MultiFab src_mf(src_ba, dmap, NUM_VARS, num_extra_ghosts,
                               amrex::MFInfo(), *m_factory);

        // Fill the multifab with the needed state data including the ghost
        // cells
        FillPatch(*this, src_mf, num_extra_ghosts, time, state_idx,
                  derive_scomp, derive_ncomp);

        const auto &src_arrays = src_mf.const_arrays();

        Potential potential(simParams().potential_params);
        ScalarFieldWithPotential scalar_field(potential);

        if (name == "constraints")
        {
            const auto &out_arrays = multifab.arrays();
            int iham               = dcomp;
            Interval imom = Interval(dcomp + 1, dcomp + AMREX_SPACEDIM);
            MatterConstraints<ScalarFieldWithPotential> constraints(
                scalar_field, Geom().CellSize(0), simParams().G_Newton, iham,
                imom);
            amrex::ParallelFor(
                multifab, multifab.nGrowVect(),
                [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k) noexcept {
                    constraints.compute(i, j, k, out_arrays[box_no],
                                        src_arrays[box_no]);
                });
        }
        else if (name == "Weyl4")
        {
            const auto &out_arrays = multifab.arrays();

            MatterWeyl4<ScalarFieldWithPotential> weyl4(
                scalar_field, simParams().extraction_params.center,
                Geom().CellSize(0), dcomp, simParams().formulation,
                simParams().G_Newton);

            amrex::ParallelFor(
                multifab, multifab.nGrowVect(),
                [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k) noexcept {
                    weyl4.compute(i, j, k, out_arrays[box_no],
                                  src_arrays[box_no]);
                });
        }
        else
        {
            amrex::Abort("Unknown derived variable");
        }
    }
    else
    {
        amrex::Abort("Unknown derived variable");
    }
    amrex::Gpu::streamSynchronize();
}

void ScalarFieldLevel::specificPostTimeStep(amrex::Real dt, int restart_time)
{
	BL_PROFILE("ScalarFieldLevel::specificPostTimeStep");

	amrex::MultiFab &state_new = get_new_data(State_Type);
	const auto cur_time        = get_state_data(State_Type).curTime();

	auto first_step = (cur_time == 0);
	const double volume = Geom().ProbSize(); // (!!) unitless volume
    const int nghost = simParams().num_ghosts;

	const double phi_avg = state_new.sum(c_phi)/volume;
	const double Pi_avg = state_new.sum(c_Pi)/volume;
    const double chi_avg = state_new.sum(c_chi)/volume;
	const double scale_fact_avg = 1./sqrt(chi_avg);
	const double Hubble_fact_avg = -state_new.sum(c_K)/volume/3.;
	const double lapse_avg = state_new.sum(c_lapse)/volume;

    const amrex::BoxArray& ba = state_new.boxArray();
    const amrex::DistributionMapping& dm = state_new.DistributionMap();
    int ncomp = state_new.nComp();
    int ngrow = state_new.nGrow();

    amrex::MultiFab phi_alias(ba, dm, ncomp, ngrow);
    amrex::MultiFab chi_alias(ba, dm, ncomp, ngrow);
    amrex::MultiFab::Copy(phi_alias, state_new, c_phi, c_phi, 1, ngrow);
    amrex::MultiFab::Copy(chi_alias, state_new, c_chi, c_chi, 1, ngrow);

    if(phi_alias.empty()) { amrex::Error("ScalarFieldLevel::specificPostTimeStep Copy failed"); }
    else if(chi_alias.empty()) { amrex::Error("ScalarFieldLevel::specificPostTimeStep Copy failed"); }

    Multiply(phi_alias, phi_alias, c_phi, c_phi, 1, nghost);
    Multiply(chi_alias, chi_alias, c_chi, c_chi, 1, nghost);
    const double phi_var = phi_alias.sum(c_phi)/volume - std::pow(phi_avg, 2.);
    const double chi_var = chi_alias.sum(c_chi)/volume - std::pow(chi_avg, 2.);

	SmallDataIO means_file(simParams().data_path+"means_file", dt, cur_time, restart_time, SmallDataIO::APPEND, first_step, ".dat");
	means_file.remove_duplicate_time_data(); // removes any duplicate data from previous run (for checkpointing)

    if(first_step) 
    {
        //constrs_file.write_header_line({"HamMean","HamSTD","HamAbsMean","HamNormMean","HamNormSTD","MomBar","MomAAD"});
        means_file.write_header_line({"PhiMean","PhiVar","PiMean","ScaleFactMean","ChiMean","ChiVar","HubbleMean","LapseMean"});
    }

    //constrs_file.write_time_data_line({hamBar, sqrt(hamVar), hamAbsBar, hamNormBar, sqrt(hamNormVar), momBar, momAAD});
    means_file.write_time_data_line({phi_avg, phi_var, Pi_avg, scale_fact_avg, chi_avg, chi_var, Hubble_fact_avg, lapse_avg});
}

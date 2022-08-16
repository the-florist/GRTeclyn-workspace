/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#include "BinaryBHLevel.hpp"
#include "AMRReductions.hpp"
#include "BinaryBH.hpp"
#include "BoxLoops.hpp"
#include "CCZ4RHS.hpp"
#include "ChiExtractionTaggingCriterion.hpp"
#include "ChiPunctureExtractionTaggingCriterion.hpp"
#include "ComputePack.hpp"
#include "NanCheck.hpp"
#include "NewConstraints.hpp"
#include "PositiveChiAndAlpha.hpp"
#include "PunctureTracker.hpp"
#include "SetValue.hpp"
#include "SixthOrderDerivatives.hpp"
#include "SmallDataIO.hpp"
#include "TraceARemoval.hpp"
#include "TwoPuncturesInitialData.hpp"
#include "Weyl4.hpp"
#include "WeylExtraction.hpp"

// Things to do during the advance step after RK4 steps
void BinaryBHLevel::specificAdvance()
{
#if 0
//xxxxx
    // Enforce the trace free A_ij condition and positive chi and alpha
    BoxLoops::loop(make_compute_pack(TraceARemoval(), PositiveChiAndAlpha()),
                   m_state_new, m_state_new, INCLUDE_GHOST_CELLS);

    // Check for nan's
    if (m_p.nan_check)
        BoxLoops::loop(NanCheck("NaNCheck in specific Advance: "), m_state_new,
                       m_state_new, EXCLUDE_GHOST_CELLS, disable_simd());
#endif
}

// This initial data uses an approximation for the metric which
// is valid for small boosts
void BinaryBHLevel::initData()
{
    BL_PROFILE("BinaryBHLevel::initialData");
    if (m_verbosity)
        amrex::Print() << "BinaryBHLevel::initialData " << Level() << std::endl;
#ifdef USE_TWOPUNCTURES
    // xxxx USE_TWOPUNCTURES todo
    TwoPuncturesInitialData two_punctures_initial_data(
        m_dx, m_p.center, m_tp_amr.m_two_punctures);
    // Can't use simd with this initial data
    BoxLoops::loop(two_punctures_initial_data, m_state_new, m_state_new,
                   INCLUDE_GHOST_CELLS, disable_simd());
#else
    // Set up the compute class for the BinaryBH initial data
    BinaryBH binary(simParams().bh1_params, simParams().bh2_params,
                    Geom().CellSize(0));

    // First set everything to zero (to avoid undefinded values in constraints)
    // then calculate initial data
    amrex::MultiFab& state = get_new_data(State_Type);
    auto const& arrs = state.arrays();
    amrex::ParallelFor(state, state.nGrowVect(),
    [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k)
    {
        amrex::Array4<amrex::Real> const& a = arrs[box_no];
        for (int n = 0; n < a.nComp(); ++n) {
            a(i,j,k,n) = 0.;
        }
        // xxxxx this is hack that's going to very slow.
        // We need to modify Bnindary BH to take Array4.
	amrex::FArrayBox fab(amrex::Box(a), a.nComp(), a.dataPtr());
        BoxPointers box_pointers(fab,fab);
        binary.compute(Cell<double>(amrex::IntVect(i,j,k), box_pointers));
    });
#endif
}

// Calculate RHS during RK4 substeps
void BinaryBHLevel::specificEvalRHS(amrex::MultiFab const & a_soln,
                                    amrex::MultiFab &a_rhs,
                                    const double a_time)
{
    // Enforce positive chi and alpha and trace free A
    

//    BoxLoops::loop(make_compute_pack(TraceARemoval(), PositiveChiAndAlpha()),
//                   a_soln, a_soln, INCLUDE_GHOST_CELLS);

    // Calculate CCZ4 right hand side
    if (simParams().max_spatial_derivative_order == 4)
    {
//        BoxLoops::loop(CCZ4RHS<MovingPunctureGauge, FourthOrderDerivatives>(
//                           m_p.ccz4_params, m_dx, m_p.sigma, m_p.formulation),
//                       a_soln, a_rhs, EXCLUDE_GHOST_CELLS);
    }
    else if (simParams().max_spatial_derivative_order == 6)
    {
//        BoxLoops::loop(CCZ4RHS<MovingPunctureGauge, SixthOrderDerivatives>(
//                           m_p.ccz4_params, m_dx, m_p.sigma, m_p.formulation),
//                       a_soln, a_rhs, EXCLUDE_GHOST_CELLS);
    }
}

// enforce trace removal during RK4 substeps
void BinaryBHLevel::specificUpdateODE(GRLevelData &a_soln,
                                      const GRLevelData &a_rhs, amrex::Real a_dt)
{
    // Enforce the trace free A_ij condition
    BoxLoops::loop(TraceARemoval(), a_soln, a_soln, INCLUDE_GHOST_CELLS);
}

void BinaryBHLevel::preTagCells()
{
    // We only use chi in the tagging criterion so only fill the ghosts for chi
//xxxxx    fillAllGhosts(VariableType::evolution, Interval(c_chi, c_chi));
}

// specify the cells to tag
void BinaryBHLevel::computeTaggingCriterion(amrex::FArrayBox &tagging_criterion,
                                            const amrex::FArrayBox &current_state)
{
#if 0
//xxxxx
    if (m_p.track_punctures)
    {
        std::vector<double> puncture_masses;
#ifdef USE_TWOPUNCTURES
        // use calculated bare masses from TwoPunctures
        puncture_masses = {m_tp_amr.m_two_punctures.mm,
                           m_tp_amr.m_two_punctures.mp};
#else
        puncture_masses = {m_p.bh1_params.mass, m_p.bh2_params.mass};
#endif /* USE_TWOPUNCTURES */
        auto puncture_coords =
            m_bh_amr.m_puncture_tracker.get_puncture_coords();
        BoxLoops::loop(ChiPunctureExtractionTaggingCriterion(
                           m_dx, m_level, m_p.max_level, m_p.extraction_params,
                           puncture_coords, m_p.activate_extraction,
                           m_p.track_punctures, puncture_masses),
                       current_state, tagging_criterion);
    }
    else
    {
        BoxLoops::loop(ChiExtractionTaggingCriterion(m_dx, m_level,
                                                     m_p.extraction_params,
                                                     m_p.activate_extraction),
                       current_state, tagging_criterion);
    }
#endif
}

void BinaryBHLevel::specificPostTimeStep()
{
#if 0
//xxxxx
    BL_PROFILE("BinaryBHLevel::specificPostTimeStep");

    bool first_step =
        (m_time == 0.); // this form is used when 'specificPostTimeStep' was
                        // called during setup at t=0 from Main
    // bool first_step = (m_time == m_dt); // if not called in Main

    if (m_p.activate_extraction == 1)
    {
        int min_level = m_p.extraction_params.min_extraction_level();
        bool calculate_weyl = at_level_timestep_multiple(min_level);
        if (calculate_weyl)
        {
            // Populate the Weyl Scalar values on the grid
            fillAllGhosts();
            BoxLoops::loop(
                Weyl4(m_p.extraction_params.center, m_dx, m_p.formulation),
                m_state_new, m_state_diagnostics, EXCLUDE_GHOST_CELLS);

            // Do the extraction on the min extraction level
            if (m_level == min_level)
            {
                BL_PROFILE("WeylExtraction");
                // Now refresh the interpolator and do the interpolation
                // fill ghosts manually to minimise communication
                bool fill_ghosts = false;
                m_gr_amr.m_interpolator->refresh(fill_ghosts);
                m_gr_amr.fill_multilevel_ghosts(
                    VariableType::diagnostic, Interval(c_Weyl4_Re, c_Weyl4_Im),
                    min_level);
                WeylExtraction my_extraction(m_p.extraction_params, m_dt,
                                             m_time, first_step,
                                             m_restart_time);
                my_extraction.execute_query(m_gr_amr.m_interpolator);
            }
        }
    }

    if (m_p.calculate_constraint_norms)
    {
        fillAllGhosts();
        BoxLoops::loop(Constraints(m_dx, c_Ham, Interval(c_Mom1, c_Mom3)),
                       m_state_new, m_state_diagnostics, EXCLUDE_GHOST_CELLS);
        if (m_level == 0)
        {
            AMRReductions<VariableType::diagnostic> amr_reductions(m_gr_amr);
            double L2_Ham = amr_reductions.norm(c_Ham);
            double L2_Mom = amr_reductions.norm(Interval(c_Mom1, c_Mom3));
            SmallDataIO constraints_file(m_p.data_path + "constraint_norms",
                                         m_dt, m_time, m_restart_time,
                                         SmallDataIO::APPEND, first_step);
            constraints_file.remove_duplicate_time_data();
            if (first_step)
            {
                constraints_file.write_header_line({"L^2_Ham", "L^2_Mom"});
            }
            constraints_file.write_time_data_line({L2_Ham, L2_Mom});
        }
    }

    // do puncture tracking on requested level
    if (m_p.track_punctures && m_level == m_p.puncture_tracking_level)
    {
        BL_PROFILE("PunctureTracking");
        // only do the write out for every coarsest level timestep
        int coarsest_level = 0;
        bool write_punctures = at_level_timestep_multiple(coarsest_level);
        m_bh_amr.m_puncture_tracker.execute_tracking(m_time, m_restart_time,
                                                     m_dt, write_punctures);
    }
#endif
}

#ifdef AMREX_USE_HDF5
// Things to do before a plot level - need to calculate the Weyl scalars
void BinaryBHLevel::prePlotLevel()
{
    fillAllGhosts();
    if (m_p.activate_extraction == 1)
    {
        BoxLoops::loop(
            make_compute_pack(
                Weyl4(m_p.extraction_params.center, m_dx, m_p.formulation),
                Constraints(m_dx, c_Ham, Interval(c_Mom1, c_Mom3))),
            m_state_new, m_state_diagnostics, EXCLUDE_GHOST_CELLS);
    }
}
#endif /* AMREX_USE_HDF5 */

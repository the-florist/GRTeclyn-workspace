/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#include "GRAMRLevel.hpp"
#include "NullBCFill.hpp"

void GRAMRLevel::stateVariableSetUp()
{
    const int nghost = simParams().num_ghosts;
    desc_lst.addDescriptor(State_Type, amrex::IndexType::TheCellType(),
                           amrex::StateDescriptor::Point, nghost, NUM_VARS,
                           &amrex::quartic_interp);

    BoundaryConditions::params_t bparms = simParams().boundary_params;
    BoundaryConditions boundary_conditions;
    boundary_conditions.define(simParams().center, bparms,
                               amrex::DefaultGeometry(), nghost);

    amrex::Vector<amrex::BCRec> bcs(NUM_VARS);
    for (int icomp = 0; icomp < NUM_VARS; ++icomp)
    {
        auto &bc = bcs[icomp];
        for (amrex::OrientationIter oit; oit.isValid(); ++oit)
        {
            amrex::Orientation face = oit();
            const int idim          = face.coordDir();
            const int bctype = boundary_conditions.get_boundary_condition(face);
            if (amrex::DefaultGeometry().isPeriodic(idim))
            {
                bc.set(face, amrex::BCType::int_dir);
            }
            else if (bctype == BoundaryConditions::STATIC_BC ||
                     bctype == BoundaryConditions::SOMMERFELD_BC ||
                     bctype == BoundaryConditions::MIXED_BC)
            {
                bc.set(face, amrex::BCType::foextrap);
            }
            else if (bctype == BoundaryConditions::REFLECTIVE_BC)
            {
                int parity =
                    BoundaryConditions::get_state_var_parity(icomp, idim);
                if (parity == 1)
                {
                    bc.set(face, amrex::BCType::reflect_even);
                }
                else
                {
                    bc.set(face, amrex::BCType::reflect_odd);
                }
            }
            else if (bctype == BoundaryConditions::EXTRAPOLATING_BC)
            {
                amrex::Abort("xxxxx EXTRAPOLATING_BC todo");
            }
            else
            {
                amrex::Abort("Unknow BC type " + std::to_string(bctype));
            }
        }
    }

    amrex::StateDescriptor::BndryFunc bndryfunc(null_bc_fill);
    bndryfunc.setRunOnGPU(true); // Run the bc function on gpu.

    desc_lst.setComponent(State_Type, 0, StateVariables::names, bcs, bndryfunc);
}

void GRAMRLevel::variableCleanUp()
{
    desc_lst.clear();
    derive_lst.clear();
}

GRAMRLevel::GRAMRLevel() = default;

GRAMRLevel::GRAMRLevel(amrex::Amr &papa, int lev, const amrex::Geometry &geom,
                       const amrex::BoxArray &box_array,
                       const amrex::DistributionMapping &distribution_mapping,
                       amrex::Real time)
    : amrex::AmrLevel(papa, lev, geom, box_array, distribution_mapping, time),
      m_num_ghosts(simParams().num_ghosts)
{

    m_boundaries.define(simParams().center, simParams().boundary_params, geom,
                        m_num_ghosts);
}

GRAMRLevel::~GRAMRLevel() = default;

const SimulationParameters &GRAMRLevel::simParams()
{
    return GRAMR::get_simulation_parameters();
}

GRAMR *GRAMRLevel::get_gramr_ptr()
{
    if (m_gramr_ptr == nullptr)
    {
        if (parent == nullptr)
        {
            amrex::Abort("AmrLevel::parent is null");
        }
        m_gramr_ptr = dynamic_cast<GRAMR *>(parent);
    }
    return m_gramr_ptr;
}

void GRAMRLevel::computeInitialDt(
    int finest_level, int /*sub_cycle*/, amrex::Vector<int> & /*n_cycle*/,
    const amrex::Vector<amrex::IntVect> & /*ref_ratio*/,
    amrex::Vector<amrex::Real> &dt_level, amrex::Real /*stop_time*/)
{
    // Level 0 will do it for all levels
    if (Level() == 0)
    {
        double dt_multiplier = simParams().dt_multiplier;
        for (int i = 0; i <= finest_level; ++i)
        {
            dt_level[i] = dt_multiplier * parent->Geom(i).CellSize(0);
        }
    }
}

void GRAMRLevel::computeNewDt(
    int finest_level, int /*sub_cycle*/, amrex::Vector<int> & /*n_cycle*/,
    const amrex::Vector<amrex::IntVect> & /*ref_ratio*/,
    amrex::Vector<amrex::Real> &dt_min, amrex::Vector<amrex::Real> &dt_level,
    amrex::Real /*stop_time*/, int /*post_regrid_flag*/)
{
    // This is called at the end of a coarse time step
    // Level 0 will do it for all levels
    if (Level() == 0)
    {
        double dt_multiplier = simParams().dt_multiplier;
        for (int i = 0; i <= finest_level; ++i)
        {
            dt_min[i] = dt_level[i] =
                dt_multiplier * parent->Geom(i).CellSize(0);
        }
    }
}

amrex::Real GRAMRLevel::advance(amrex::Real time, amrex::Real dt, int iteration,
                                int ncycle)
{
    BL_PROFILE("GRAMRLevel::advance()");
    double seconds_per_hour = 3600;
    double evolution_speed  = (time - get_gramr_ptr()->get_restart_time()) *
                             seconds_per_hour /
                             get_gramr_ptr()->get_walltime_since_start();
    amrex::Print() << "[Level " << Level() << " step "
                   << parent->levelSteps(Level()) + 1
                   << "] average evolution speed = " << evolution_speed
                   << " code units/h\n";

    for (int k = 0; k < NUM_STATE_TYPE; k++)
    {
        state[k].allocOldData();
        state[k].swapTimeLevels(dt);
    }

    amrex::AmrLevel::RK(
        4, State_Type, time, dt, iteration, ncycle,
        [&](int /*stage*/, amrex::MultiFab &rhs, const amrex::MultiFab &soln,
            amrex::Real t, amrex::Real /*dtsub*/)
        {
            // NOLINTNEXTLINE(cppcoreguidelines-pro-type-const-cast)
            specificEvalRHS(const_cast<amrex::MultiFab &>(soln), rhs, t);
            m_boundaries.apply_sommerfeld_boundaries(rhs, soln);

            amrex::Gpu::streamSynchronize();
        },
        [&](int /*stage*/, amrex::MultiFab &soln) { specificUpdateODE(soln); });

    specificAdvance();

    return dt;
}

void GRAMRLevel::post_timestep(int /*iteration*/)
{
    BL_PROFILE("GRAMRLevel::post_timestep()");
    const int lev = Level();
    if (lev < parent->finestLevel())
    {
        auto &fine_level        = parent->getLevel(Level() + 1);
        amrex::MultiFab &S_fine = fine_level.get_new_data(State_Type);
        amrex::MultiFab &S_crse = this->get_new_data(State_Type);
        amrex::Real t           = get_state_data(State_Type).curTime();

        amrex::IntVect ratio = parent->refRatio(lev);
        AMREX_ASSERT(ratio == 2 || ratio == 4);
        if (ratio == 2)
        {
            // Need to fill one ghost cell for the high-order interpolation
            // below
            FillPatch(fine_level, S_fine, 1, t, State_Type, 0, S_fine.nComp());
        }

        FourthOrderInterpFromFineToCoarse(S_crse, 0, NUM_VARS, S_fine, ratio);
    }

    amrex::Real dt = parent->dtLevel(level);
    int restart_time = get_gramr_ptr()->get_restart_time();
    specificPostTimeStep(dt, restart_time);
}

void GRAMRLevel::post_regrid(int /*lbase*/, int /*new_finest*/)
{
    // xxxxx Do we need to do anything after regrid?
}

void GRAMRLevel::post_init(amrex::Real /*stop_time*/)
{
    if (Level() == 0)
    {
        get_gramr_ptr()->set_restart_time(get_gramr_ptr()->cumTime());
    }

    amrex::Real dt = parent->dtLevel(level);
    int restart_time = get_gramr_ptr()->get_restart_time();
    specificPostTimeStep(dt, restart_time);
}

void GRAMRLevel::post_restart()
{
    if (Level() == 0)
    {
        get_gramr_ptr()->set_restart_time(get_gramr_ptr()->cumTime());
    }
}

void GRAMRLevel::init(amrex::AmrLevel &old)
{
    BL_PROFILE("GRAMRLevel::init()");
    amrex::Real dt_new    = parent->dtLevel(level);
    amrex::Real cur_time  = old.get_state_data(State_Type).curTime();
    amrex::Real prev_time = old.get_state_data(State_Type).prevTime();
    amrex::Real dt_old    = cur_time - prev_time;
    setTimeLevel(cur_time, dt_old, dt_new);

    amrex::MultiFab &S_new = get_new_data(State_Type);
    FillPatch(old, S_new, 0, cur_time, State_Type, 0, S_new.nComp());
}

void GRAMRLevel::init()
{
    BL_PROFILE("GRAMRLevel::init()");
    amrex::Real dt = parent->dtLevel(level);
    const auto &coarse_state =
        parent->getLevel(level - 1).get_state_data(State_Type);
    amrex::Real cur_time  = coarse_state.curTime();
    amrex::Real prev_time = coarse_state.prevTime();
    amrex::Real dt_old =
        (cur_time - prev_time) /
        static_cast<amrex::Real>(parent->MaxRefRatio(level - 1));
    setTimeLevel(cur_time, dt_old, dt);

    amrex::MultiFab &S_new = get_new_data(State_Type);
    FillCoarsePatch(S_new, 0, cur_time, State_Type, 0, S_new.nComp());
}

void GRAMRLevel::writePlotFilePre(const std::string & /*dir*/,
                                  std::ostream & /*os*/)
{
    m_is_writing_plotfile = true;
    // auto &state_new       = get_new_data(State_Type);
    // FillPatch(*this, state_new, state_new.nGrow(),
    //           get_state_data(State_Type).curTime(), State_Type, 0,
    //           state_new.nComp());
}

void GRAMRLevel::writePlotFilePost(const std::string & /*dir*/,
                                   std::ostream & /*os*/)
{
    m_is_writing_plotfile = false;
}

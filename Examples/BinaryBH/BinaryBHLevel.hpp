/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef BINARYBHLEVEL_HPP_
#define BINARYBHLEVEL_HPP_

#include "DefaultLevelFactory.hpp"
#include "GRAMRLevel.hpp"
// TPAMR.hpp includes BHAMR.hpp
#include "TPAMR.hpp"

class BinaryBHLevel : public GRAMRLevel
{
  public:
    static void variableSetUp();

    // Inherit the contructors from GRAMRLevel
    using GRAMRLevel::GRAMRLevel;

    /// Things to do at every full timestep
    ///(might include several substeps, e.g. in RK4)
    void specificAdvance() override;

    /// Initial data calculation
    void initData() override;

    /// Calculation of the right hand side for the time stepping
    void specificEvalRHS(amrex::MultiFab &a_soln, amrex::MultiFab &a_rhs,
                         const double a_time) override;

    /// Things to do after dt*rhs has been added to the solution
    void specificUpdateODE(amrex::MultiFab &a_soln) override;

    // to do post each time step on every level
    void specificPostTimeStep() override;

    void errorEst(amrex::TagBoxArray &tag_box_array, int clearval, int tagval,
                  amrex::Real time, int n_error_buf = 0, int ngrow = 0) final;
};

#endif /* BINARYBHLEVEL_HPP_ */

/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef SCALARFIELDLEVEL_HPP_
#define SCALARFIELDLEVEL_HPP_

#include "DefaultLevelBld.hpp"
#include "FourthOrderDerivatives.hpp"
#include "GRAmrLevel.hpp"
#include "Potential.hpp"
#include "ScalarField.hpp"
#include "SixthOrderDerivatives.hpp"

/// Evolution level for a real scalar field minimally coupled to gravity.
class ScalarFieldLevel : public GRAmrLevel
{
  public:

    // Inherit the contructors from GRAmrLevel
    using GRAmrLevel::GRAmrLevel;

    template <class deriv_t = FourthOrderDerivatives>
    using ScalarFieldWithPotential = ScalarField<Potential, deriv_t>;

    static void variableSetUp();

    //! Things to do at the end of the advance step, after RK4 calculation
    void specific_advance() override;

    //! Initialize data for the field and metric variables
    void initData() override;

    //! RHS routines used at each RK4 step
    void specific_eval_rhs(amrex::MultiFab &a_soln, amrex::MultiFab &a_rhs,
                         const amrex::Real a_time) override;

    //! Things to do in UpdateODE step, after soln + rhs update
    void specific_update_ode(amrex::MultiFab &a_soln) override;

    void specific_post_timestep() override;

    //! Run the post-timestep processing once on the initial data too
    void specific_post_init() override;

    void tag_cells(amrex::TagBoxArray &a_tag_box_array,
                amrex::Real a_regrid_threshold) final;
};

#endif /* SCALARFIELDLEVEL_HPP_ */

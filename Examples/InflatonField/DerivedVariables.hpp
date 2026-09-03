/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef DERIVEDVARIABLES_HPP_
#define DERIVEDVARIABLES_HPP_

#include "FilesystemTools.hpp"
#include "InflatonUtils.hpp"
#include "SmallDataIO.hpp"

#include <AMReX_FFT.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_Print.H>
#include <AMReX_Vector.H>

class DerivedVariables
{
  public:
    //! Name of the derive record. This, not the individual variable names
    //! below, is what amr.derive_plot_vars is matched against.
    static inline const std::string name = "InflationFields";

    // Names of diagnostic variables
    static inline const amrex::Vector<std::string> var_names{"R", "hplus",
                                                             "hcross"};

    static void set_up(int a_state_index);

    // Derive callback (amrex::DeriveFuncMF) that fills plotfile output
    static void compute_mf(amrex::MultiFab &out_mf, int dcomp, int ncomp,
                           const amrex::MultiFab &src_mf,
                           const amrex::Geometry &geomdata, amrex::Real time,
                           const int *bcrec, int level);

  private:
    InflatonUtils m_utils;

    [[nodiscard]] const InflatonParameters &params() const
    {
        return m_utils.m_params;
    }

    void extract_hs_and_R(amrex::MultiFab &hs, amrex::MultiFab &R,
                          const amrex::MultiFab &state);
};

#include "DerivedVariables.impl.hpp"

#endif /* DERIVEDVARIABLES_HPP_ */
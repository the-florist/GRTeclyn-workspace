/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef INFLATIONEXTRACTION_HPP_
#define INFLATIONEXTRACTION_HPP_

#include "Config.hpp"
#include "FilesystemTools.hpp"
#include "SmallDataIO.hpp"

#include <AMReX_FFT.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_Print.H>
#include <AMReX_Vector.H>

class InflationExtraction
{
  public:
    // Names of diagnostic variables
    static inline const amrex::Vector<std::string> var_names{"R", "hplus",
                                                             "hcross"};

    // Constructor used in extraction of diagnostics
    InflationExtraction() { m_inflaton_methods.fill_params(); }

    static void set_up(int a_state_index);

    // Derive callback (amrex::DeriveFuncMF) that fills plotfile output
    static void compute_mf(amrex::MultiFab &out_mf, int dcomp, int ncomp,
                           const amrex::MultiFab &src_mf,
                           const amrex::Geometry &geomdata, amrex::Real time,
                           const int *bcrec, int level);

  private:
    Config m_inflaton_methods;

    void extract_hs_and_R(amrex::MultiFab &hs, amrex::MultiFab &R,
                          const amrex::MultiFab &state);
};

#include "InflationExtraction.impl.hpp"

#endif /* INFLATIONEXTRACTION_HPP_ */
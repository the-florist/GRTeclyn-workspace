/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef INFLATONFIELDINIT_HPP_
#define INFLATONFIELDINIT_HPP_

#include "GRParmParse.hpp"
#include "InflationConfig.hpp"
#include "StateVariables.hpp"
#include "TensorTests.hpp"

#include <AMReX_FFT.H>
#include <AMReX_GpuContainers.H>

class InflatonFieldInit
{
  protected:
    InflationConfig inflt_methods;

  public:
    // Constructor used when initialising stochastic fields
    InflatonFieldInit()
    {
        // Read the parameters once from the global GRParmParse table
        inflt_methods.fill_params();
    }

    void init(amrex::MultiFab &state);

  private:
    enum class FieldType
    {
        Scalar,
        Tensor
    };
    enum class WhichField
    {
        Amplitude = 0,
        Velocity  = 1
    };
    enum class BSSNFields
    {
        Phi = 0,
        Pi  = 1,
        Chi = 2,
        K   = 3
    };

    AMREX_GPU_HOST_DEVICE static amrex::GpuComplex<amrex::Real>
    calculate_mode_function(const InflationConfig &cfg, const amrex::Real km,
                            const FieldType field_type,
                            const WhichField which_field);

    AMREX_GPU_HOST_DEVICE static amrex::GpuComplex<amrex::Real>
    calculate_random_field(const InflationConfig &cfg, const amrex::IntVect iv,
                           const amrex::Real rand_amp,
                           const amrex::Real rand_phase,
                           const FieldType field_type,
                           const WhichField which_field);

    void generate_fourier_realisation(amrex::cMultiFab &hij_k,
                                      amrex::cMultiFab &Aij_k,
                                      amrex::cMultiFab &scalar_fields_k);

    void add_perturbations_to_state(amrex::MultiFab &state,
                                    amrex::MultiFab &hij_x,
                                    amrex::MultiFab &Aij_x,
                                    amrex::MultiFab &scalar_fields_x,
                                    const int dN);
};

#include "InflatonFieldInit.impl.hpp"

#endif /* INFLATONFIELDINIT_HPP_ */

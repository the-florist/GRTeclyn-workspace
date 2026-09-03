/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef SCALARTENSORINIT_HPP_
#define SCALARTENSORINIT_HPP_

#include "GRParmParse.hpp"
#include "InflatonUtils.hpp"
#include "StateVariables.hpp"

#include <AMReX_FFT.H>
#include <AMReX_GpuContainers.H>

class ScalarTensorInit
{
  protected:
    InflatonUtils m_utils;

    const InflatonParameters &params() const { return m_utils.m_params; }

  public:
    // Constructor used when initialising stochastic fields
    ScalarTensorInit() { ; }

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
    calculate_mode_function(const InflatonParameters &d_params,
                            const amrex::Real km, const FieldType field_type,
                            const WhichField which_field);

    AMREX_GPU_HOST_DEVICE static amrex::GpuComplex<amrex::Real>
    calculate_random_field(const InflatonUtils &cfg,
                           const InflatonParameters &d_params,
                           const amrex::IntVect iv, const amrex::Real rand_amp,
                           const amrex::Real rand_phase,
                           const FieldType field_type,
                           const WhichField which_field);

    AMREX_GPU_HOST_DEVICE inline void 
    ScalarTensorInit::convert_R_to_BSSN_scalars(const InflatonUtils &cfg,
                                                const InflatonParameters &params,
                                                amrex::cMultiFab &R_and_dR,
                                                amrex::cMultiFab &bssn_scalars)

    void generate_fourier_realisation(amrex::cMultiFab &hij_k,
                                      amrex::cMultiFab &Aij_k,
                                      amrex::cMultiFab &scalar_fields_k);

    void add_perturbations_to_state(amrex::MultiFab &state,
                                    amrex::MultiFab &hij_x,
                                    amrex::MultiFab &Aij_x,
                                    amrex::MultiFab &scalar_fields_x,
                                    const int dN);
};

#include "ScalarTensorInit.impl.hpp"

#endif /* SCALARTENSORINIT_HPP_ */

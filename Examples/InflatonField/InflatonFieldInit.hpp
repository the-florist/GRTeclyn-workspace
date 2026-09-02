/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef RANDOMFIELDINIT_HPP_
#define RANDOMFIELDINIT_HPP_

#include "GRParmParse.hpp"
#include "InflationConfig.hpp"
#include "StateVariables.hpp"
#include "TensorTests.hpp"

#include <AMReX_FFT.H>
#include <AMReX_GpuContainers.H>

class RandomFieldInit
{
  protected:
    InflationConfig inflt_methods;

  public:
    // Constructor used when initialising stochastic fields
    RandomFieldInit()
    {
        // Read the parameters once from the global GRParmParse table
        inflt_methods.fill_params();

        // Flatten and upload the STOIIC spectra to the device
        upload_stoiic_to_device();
    }

    void init(amrex::MultiFab &state);

  private:
    enum class FieldType
    {
        Scalar,
        Tensor
    };
    enum class ScalarField
    {
        Phi = 0,
        Pi  = 1,
        Chi = 2,
        K   = 3
    };
    enum class TensorField
    {
        Amplitude = 0,
        Velocity  = 1
    };

    // Device storage backing the STOIIC pointers in inflt_methods
    amrex::Gpu::DeviceVector<amrex::Real> m_init_k_d;
    amrex::Gpu::DeviceVector<amrex::Real> m_scalar_ps_d;
    amrex::Gpu::DeviceVector<amrex::Real> m_tensor_ps_d;

    void upload_stoiic_to_device();

    AMREX_GPU_HOST_DEVICE static amrex::GpuComplex<amrex::Real>
    calculate_mode_function(const InflationConfig &cfg, const amrex::Real km,
                            const TensorField field_selector);

    AMREX_GPU_HOST_DEVICE static amrex::GpuComplex<amrex::Real>
    find_in_stoiic(const InflationConfig &cfg, const amrex::Real km,
                   const FieldType field_type, const auto field_selector);

    AMREX_GPU_HOST_DEVICE static amrex::GpuComplex<amrex::Real>
    calculate_random_field(const InflationConfig &cfg, const amrex::IntVect iv,
                           const amrex::Real rand_amp,
                           const amrex::Real rand_phase,
                           const FieldType field_type,
                           const auto field_selector);

    amrex::Real find_precision_loss(amrex::MultiFab &field, const int comp,
                                    const amrex::Real bkgd);

    void generate_fourier_realisation(amrex::cMultiFab &hij_k,
                                      amrex::cMultiFab &Aij_k,
                                      amrex::cMultiFab &scalar_fields_k);

    void add_perturbations_to_state(amrex::MultiFab &state,
                                    amrex::MultiFab &hij_x,
                                    amrex::MultiFab &Aij_x,
                                    amrex::MultiFab &scalar_fields_x,
                                    const int dN);
};

#include "RandomFieldInit.impl.hpp"

#endif /* RANDOMFIELDINIT_HPP_ */

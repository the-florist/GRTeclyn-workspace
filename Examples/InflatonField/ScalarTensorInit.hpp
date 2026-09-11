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

#include <vector>

// GPU-callable linear interpolator over a tabulated complex spectrum
// (modulus/phase interpolation, matching ISTORIZ's
// ComplexLinearInterpolator). Only holds raw device pointers, so it is
// trivially copyable and safe to capture by value into a device lambda; the
// amrex::Gpu::DeviceVector storage it points into must outlive its use (see
// ScalarTensorInit::generate_fourier_realisation).
struct ComplexLinearInterpolator
{
    const double *k  = nullptr;
    const double *re = nullptr;
    const double *im = nullptr;
    int n            = 0;

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::GpuComplex<amrex::Real>
    operator()(amrex::Real x) const;
};

// Copies k/re/im into device-visible storage owned by the caller (d_k/d_re/
// d_im, which must stay alive for as long as the returned interpolator is
// used) and returns a GPU-callable functor.
inline ComplexLinearInterpolator make_complex_linear_interpolator(
    const std::vector<double> &k, const std::vector<double> &re,
    const std::vector<double> &im, amrex::Gpu::DeviceVector<double> &d_k,
    amrex::Gpu::DeviceVector<double> &d_re,
    amrex::Gpu::DeviceVector<double> &d_im);

class ScalarTensorInit
{
  protected:
    InflatonUtils m_utils;

    [[nodiscard]] const InflatonParameters &params() const
    {
        return m_utils.m_params;
    }

    // Host-only STOIIC_GR/ISTORIZ spectrum table (init.use_stoiic_spectra).
    // Kept out of InflatonParameters/InflatonUtils, which are captured by
    // value into device lambdas and must stay POD.
    std::vector<double> m_spectra_k;
    std::vector<double> m_spectra_re_R;
    std::vector<double> m_spectra_im_R;
    std::vector<double> m_spectra_re_dR;
    std::vector<double> m_spectra_im_dR;

    void load_stoiic_spectra();

  public:
    // Constructor used when initialising stochastic fields
    ScalarTensorInit() { ; }

    void init(amrex::MultiFab &state);

    // nvcc requires the enclosing function of an extended __device__ lambda
    // to have public access, so the kernel-launching functions below must be
    // public even though they are implementation details of init().
    void convert_R_to_BSSN_scalars(const InflatonUtils &cfg,
                                   const InflatonParameters &d_params,
                                   const amrex::cMultiFab &R_and_dR,
                                   amrex::cMultiFab &bssn_scalars);

    void generate_fourier_realisation(amrex::cMultiFab &hij_k,
                                      amrex::cMultiFab &Aij_k,
                                      amrex::cMultiFab &scalar_fields_k);

    void add_perturbations_to_state(amrex::MultiFab &state,
                                    amrex::MultiFab &hij_x,
                                    amrex::MultiFab &Aij_x,
                                    amrex::MultiFab &scalar_fields_x,
                                    const int dn_ratio);

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
                            const amrex::Real kmag, const FieldType field_type,
                            const WhichField which_field);

    AMREX_GPU_HOST_DEVICE static amrex::GpuComplex<amrex::Real>
    calculate_random_field(const InflatonUtils &cfg,
                           const InflatonParameters &d_params,
                           const amrex::IntVect ivec,
                           const amrex::Real rand_amp,
                           const amrex::Real rand_phase,
                           const FieldType field_type,
                           const WhichField which_field,
                           const ComplexLinearInterpolator &interp_R,
                           const ComplexLinearInterpolator &interp_dR);
};

#include "ScalarTensorInit.impl.hpp"

#endif /* SCALARTENSORINIT_HPP_ */

/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef RANDOMFIELDINIT_HPP_
#define RANDOMFIELDINIT_HPP_

#include "InflationConfig.hpp"
#include "InitialBackgroundData.hpp"
#include "Potential.hpp"
#include "TensorTests.hpp"

#include <AMReX_FFT.H>
#include <AMReX_GpuContainers.H>

class RandomFieldInit
{
    protected:
        amrex::Real H0;
        amrex::Real phi0;

        static amrex::Real
        calc_H0(amrex::Real G, amrex::Real Pi, amrex::Real V)
        {
            return sqrt((8. * M_PI * G / 3.)*(0.5*pow(Pi, 2.) + V));
        }

    public:
                // Constructor used when initialising stochastic fields
        RandomFieldInit(const InflationConfig a_config,
                        const InitialBackgroundData::params_t bkgd_params,
                        const Potential::params_t potential_params)
                        : m_params(a_config)
        {
            // Compute background potential
            amrex::Real V, dV;
            Potential potential(potential_params);
            switch (potential_params.type)
            {
                case 1:
                    potential.quadratic(V, dV, bkgd_params.phi0);
                    break;
                case 4:
                    potential.quadratic_bump(V, dV, bkgd_params.phi0);
                    break;
                case 8:
                    potential.USR(V, dV, bkgd_params.phi0);
                    break;
                case 9:
                    potential.monodromy(V, dV, bkgd_params.phi0);
                    break;
                case 10:
                    potential.punctuated(V, dV, bkgd_params.phi0);
                    break;
                case 11:
                    potential.inverted_quadratic_bump(V, dV, bkgd_params.phi0);
                    break;
                case 12:
                    potential.quadratic_step(V, dV, bkgd_params.phi0);
                    break;
                default:
                    amrex::Error("RandomFieldInit::RandomFieldInit,"
                                 "potential type not provided");
            }

            // Compute initial Hubble parameter
            H0 = calc_H0(bkgd_params.G_Newton, bkgd_params.Pi0, V);
            phi0 = bkgd_params.phi0;

            // Flatten and upload the STOIIC spectra to the device
            upload_stoiic_to_device();
        }

        void init(amrex::MultiFab &state);

    private:
        enum class FieldType   { Scalar, Tensor };
        enum class ScalarField { Phi = 0, Pi = 1, Chi = 2, K = 3 };
        enum class TensorField { Amplitude = 0, Velocity = 1 };

        InflationConfig m_params;

        // Device storage backing the STOIIC pointers in m_params
        amrex::Gpu::DeviceVector<amrex::Real> m_init_k_d;
        amrex::Gpu::DeviceVector<amrex::Real> m_scalar_ps_d;
        amrex::Gpu::DeviceVector<amrex::Real> m_tensor_ps_d;

        void upload_stoiic_to_device();

        AMREX_GPU_HOST_DEVICE static amrex::GpuComplex<amrex::Real>
        calculate_mode_function(const InflationParams &cfg, const amrex::Real H0,
                                const amrex::Real km,
                                const TensorField field_selector);

        AMREX_GPU_HOST_DEVICE static amrex::GpuComplex<amrex::Real>
        find_in_stoiic(const InflationParams &cfg, const amrex::Real km,
                       const FieldType field_type,
                       const auto field_selector);

        AMREX_GPU_HOST_DEVICE static amrex::GpuComplex<amrex::Real>
        calculate_random_field(const InflationParams &cfg, const amrex::Real H0,
                               const amrex::IntVect iv,
                               const amrex::Real rand_amp,
                               const amrex::Real rand_phase,
                               const FieldType field_type,
                               const auto field_selector);

        amrex::Real
        find_precision_loss(amrex::MultiFab &field,
                            const int comp,
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
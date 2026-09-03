/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#if !defined(SCALARTENSORINIT_HPP_)
#error "This file should only be included via ScalarTensorInit.hpp"
#endif

#ifndef SCALARTENSORINIT_IMPL_HPP_
#define SCALARTENSORINIT_IMPL_HPP_

// Returns analytic power spectrum in modulus/argument form
AMREX_GPU_HOST_DEVICE inline amrex::GpuComplex<amrex::Real>
ScalarTensorInit::calculate_mode_function(const InflatonParameters &d_params,
                                          const amrex::Real kmag,
                                          const FieldType field_type,
                                          const WhichField which_field)
{
    // Deals with k=0 case, which is undefined if m=0
    if (kmag < InflatonUtils::tolerance)
    {
        return amrex::GpuComplex<amrex::Real>{0., 0.};
    }

    // Stores modulus and argument
    amrex::Real ms_mag = 0.;
    amrex::Real ms_arg = 0.;

    amrex::Real kpr = kmag / d_params.H0;
    if (which_field == WhichField::Amplitude) // Position mode funcion
    {
        ms_mag =
            sqrt((1.0 / kmag + d_params.H0 * d_params.H0 / pow(kmag, 3.)) /
                2. / pow(d_params.Mp, 2.));
        ms_arg =
            atan2((cos(kpr) + kpr * sin(kpr)), (kpr * cos(kpr) - sin(kpr)));
    }
    else // Velocity mode funcion
    {
        ms_mag = sqrt(kmag / 2. / pow(d_params.Mp, 2.));
        ms_arg = -atan2(cos(kpr), sin(kpr));
    }

    // Construct the mode function and return it
    amrex::GpuComplex<amrex::Real> mode_val(ms_mag * cos(ms_arg),
                                            ms_mag * sin(ms_arg));

    if (field_type == FieldType::Scalar)
    {
        mode_val /= std::sqrt(2. * d_params.epsilon_1);
    }

    if (which_field == WhichField::Amplitude)
    {
        mode_val /= d_params.init_a;
    }
    else
    {
        mode_val /= std::pow(d_params.init_a, 2.);
    }

    return mode_val;
}

// Turns analytic PS into GRF and applies window function if requested
AMREX_GPU_HOST_DEVICE inline amrex::GpuComplex<amrex::Real>
ScalarTensorInit::calculate_random_field(const InflatonUtils &cfg,
                                         const InflatonParameters &d_params,
                                         const amrex::IntVect ivec,
                                         const amrex::Real rand_amp,
                                         const amrex::Real rand_phase,
                                         const FieldType field_type,
                                         const WhichField which_field)
{
    amrex::GpuComplex<amrex::Real> value(0., 0.);
    amrex::Real kmag = cfg.get_kmag(ivec);
    value = calculate_mode_function(d_params, kmag, field_type, which_field);

    // Add stochastic perturbations
    // Make one random draw for the amplitude and phase individually
    amrex::Real rand_mod = sqrt(-2. * log(rand_amp));
    amrex::Real rand_arg = 2. * amrex::Math::pi<amrex::Real>() * rand_phase;

    // Multiply amplitude by Rayleigh draw
    value *= rand_mod;

    // Apply the random phase (assuming MS phase is accounted for)
    amrex::Real new_real =
        value.real() * cos(rand_arg) - value.imag() * sin(rand_arg);
    amrex::Real new_imag =
        value.real() * sin(rand_arg) + value.imag() * cos(rand_arg);
    amrex::GpuComplex<amrex::Real> new_value(new_real, new_imag);

    value = new_value;

    // Apply a window function if requested
    if (d_params.use_window == 1)
    {
        value *= cfg.calculate_window_function(kmag);
    }

    return value;
}

inline void ScalarTensorInit::convert_R_to_BSSN_scalars(
    const InflatonUtils &cfg, const InflatonParameters &d_params,
    const amrex::cMultiFab &R_and_dR, amrex::cMultiFab &bssn_scalars)
{
    // Refer to https://arxiv.org/abs/2502.06783 for the derivation of these
    // initial conditions. R_and_dR is indexed by WhichField (0: R, 1: dR)
    // plus a third component (2) holding dR/k^2. bssn_scalars is indexed by
    // the BSSNFields enum (0: Phi, 1: Pi, 2: Chi, 3: K).

    const amrex::Real H0        = d_params.H0;
    const amrex::Real epsilon_1 = d_params.epsilon_1;
    const amrex::Real epsilon_2 = d_params.epsilon_2;
    const amrex::Real Mp        = d_params.Mp;
    const amrex::Real init_a    = d_params.init_a;

    const amrex::Real dlnGamma = H0 * epsilon_2;

    // Phi coefficients
    const amrex::Real factor_R1 = Mp * std::sqrt(2.0 * epsilon_1);
    const amrex::Real factor_dR1invLap =
        factor_R1 * epsilon_1 * H0 * std::pow(init_a, 2.);

    // Pi coefficients. half_pi_dR_coeff is the pre-doubled value used inside
    // factor_dR2invLap's own formula; pi_dR_coeff is the (doubled, negated)
    // value that actually multiplies dR.
    const amrex::Real factor_R2        = -Mp * std::sqrt(epsilon_1 / 2.0) *
                                         (2.0 * epsilon_1 - dlnGamma / H0) * H0;
    const amrex::Real half_pi_dR_coeff = -Mp * std::sqrt(epsilon_1 / 2.0);
    const amrex::Real factor_dR2invLap = half_pi_dR_coeff * std::pow(H0, 2.) *
                                         std::pow(init_a, 2.) * epsilon_1 *
                                         (2.0 * epsilon_1 - dlnGamma / H0);
    const amrex::Real pi_dR_coeff      = -2.0 * half_pi_dR_coeff;

    // Chi coefficients
    const amrex::Real factor_dR3invLap = -2.0 * H0 * epsilon_1;

    // K coefficients
    const amrex::Real factor_R4 = 3.0 * H0 * epsilon_1;
    const amrex::Real factor_dR4invLap =
        3.0 * init_a * init_a * H0 * H0 * epsilon_1 * epsilon_1;

    constexpr int r_comp  = static_cast<int>(WhichField::Amplitude);
    constexpr int dr_comp = static_cast<int>(WhichField::Velocity);

    const auto &r_dr_arrs = R_and_dR.const_arrays();
    const auto &bssn_arrs = bssn_scalars.arrays();
    amrex::ParallelFor(
        R_and_dR,
        [=] AMREX_GPU_DEVICE(int bx, int i, int j, int k)
        {
            const amrex::IntVect iv{i, j, k};
            if (iv == amrex::IntVect{0, 0, 0})
            {
                return;
            } // Skip the zero mode

            const auto r_val       = r_dr_arrs[bx](i, j, k, r_comp);
            const auto dr_val      = r_dr_arrs[bx](i, j, k, dr_comp);
            const auto dr_over_k2 = dr_val / std::pow(cfg.get_kmag(iv), 2.);

            bssn_arrs[bx](i, j, k, static_cast<int>(BSSNFields::Phi)) =
                factor_R1 * r_val + factor_dR1invLap * dr_over_k2;
            bssn_arrs[bx](i, j, k, static_cast<int>(BSSNFields::Pi)) =
                factor_R2 * r_val + pi_dR_coeff * dr_val +
                factor_dR2invLap * dr_over_k2;
            bssn_arrs[bx](i, j, k, static_cast<int>(BSSNFields::Chi)) =
                factor_dR3invLap * dr_over_k2;
            bssn_arrs[bx](i, j, k, static_cast<int>(BSSNFields::K)) =
                factor_R4 * r_val + factor_dR4invLap * dr_over_k2;
        });
    amrex::Gpu::streamSynchronize();
}

inline void ScalarTensorInit::generate_fourier_realisation(
    amrex::cMultiFab &hij_k, amrex::cMultiFab &Aij_k,
    amrex::cMultiFab &scalar_fields_k)
{
    // Declare array to hold R and dR fields
    amrex::cMultiFab R_dR_k(scalar_fields_k.boxArray(),
                            scalar_fields_k.DistributionMap(), 2,
                            scalar_fields_k.nGrowVect(), amrex::MFInfo(),
                            scalar_fields_k.Factory());
    R_dR_k.setVal(0.0);

    // Extract arrays before ParallelFor
    const auto &hij_k_arrs        = hij_k.arrays();
    const auto &Aij_k_arrs        = Aij_k.arrays();
    const auto &scalar_field_arrs = scalar_fields_k.arrays();
    const auto &R_dR_k_arrs       = R_dR_k.arrays();

    amrex::Print() << "Starting initial condition generation/read in...\n";

    // Local copy so the kernel captures config by value, not via the host
    // `this` pointer
    const InflatonUtils cfg           = m_utils;
    const InflatonParameters d_params = params();

    amrex::ParallelFor(
        hij_k,
        [=] AMREX_GPU_DEVICE(int bx, int i, int j, int k)
        {
            amrex::IntVect iv = {i, j, k};
            amrex::GpuArray<amrex::GpuComplex<amrex::Real>, 2>
                h_mode_function; // Amp mode fn
            amrex::GpuArray<amrex::GpuComplex<amrex::Real>, 2>
                A_mode_function; // Vel mode fn

            // Uniform random draw, MPI/OpenMP safe
            const uint64_t cell_key =
                uint64_t(d_params.random_seed) *
                (645950ULL * uint64_t(iv[0]) +
                 520666ULL * uint64_t(cfg.invert_index_with_sign(iv[1])) +
                 767051ULL * uint64_t(cfg.invert_index_with_sign(iv[2])));

            // Initialise scalar sector (one random draw)
            if (d_params.scalar_init)
            {
                amrex::Real draw1 = InflatonUtils::to_unit_open(
                    InflatonUtils::splitmix64(cell_key + 0ULL));
                amrex::Real draw2 = InflatonUtils::to_unit_open(
                    InflatonUtils::splitmix64(cell_key + 1ULL));

                R_dR_k_arrs[bx](i, j, k,
                                static_cast<int>(WhichField::Amplitude)) =
                    calculate_random_field(cfg, d_params, iv, draw1, draw2,
                                           FieldType::Scalar,
                                           WhichField::Amplitude);

                R_dR_k_arrs[bx](i, j, k,
                                static_cast<int>(WhichField::Velocity)) =
                    calculate_random_field(cfg, d_params, iv, draw1, draw2,
                                           FieldType::Scalar,
                                           WhichField::Velocity);
            }

            // Initialise tensor sector (two random draws)
            if (d_params.tensor_init)
            {
                // Find the mode function realisation
                for (int p = 0; p < 2; p++)
                {
                    // One draw per polarisation field
                    amrex::Real draw1 =
                        InflatonUtils::to_unit_open(InflatonUtils::splitmix64(
                            cell_key + uint64_t(2 + 2 * p)));
                    amrex::Real draw2 =
                        InflatonUtils::to_unit_open(InflatonUtils::splitmix64(
                            cell_key + uint64_t(3 + 2 * p)));

                    h_mode_function[p] = calculate_random_field(
                        cfg, d_params, iv, draw1, draw2, FieldType::Tensor,
                        WhichField::Amplitude);

                    A_mode_function[p] = calculate_random_field(
                        cfg, d_params, iv, draw1, draw2, FieldType::Tensor,
                        WhichField::Velocity);
                }

                // Construct polarisation tensors from basis vectors
                const auto [eplus, ecross] =
                    cfg.calculate_polarisation_tensors(iv);

                // Construct Fourier space tensors
                for (int l = 0; l < 3; l++)
                    for (int p = 0; p < 3; p++)
                    {
                        hij_k_arrs[bx](i, j, k,
                                       InflatonUtils::look_up_table[l][p]) =
                            (h_mode_function[0] * eplus(l, p) +
                             h_mode_function[1] * ecross(l, p));
                        Aij_k_arrs[bx](i, j, k,
                                       InflatonUtils::look_up_table[l][p]) =
                            (A_mode_function[0] * eplus(l, p) +
                             A_mode_function[1] * ecross(l, p));
                    }
            }
        });

    amrex::Gpu::streamSynchronize();

    convert_R_to_BSSN_scalars(m_utils, params(), R_dR_k, scalar_fields_k);

    // Apply the DC and Nyquist symmetry conditions
    m_utils.apply_nyquist_conditions(hij_k);
    m_utils.apply_nyquist_conditions(Aij_k);
    m_utils.apply_nyquist_conditions(scalar_fields_k);
}

inline void ScalarTensorInit::add_perturbations_to_state(
    amrex::MultiFab &state, amrex::MultiFab &hij_x, amrex::MultiFab &Aij_x,
    amrex::MultiFab &scalar_fields_x, const int dn_ratio)
{
    static_assert(c_h33 == c_h11 + 5 && c_A33 == c_A11 + 5,
                  "add_perturbations_to_state assumes c_h11..c_h33 and "
                  "c_A11..c_A33 are consecutive metric/A_ij components");

    // Apply normalisation into physical units
    hij_x.mult(m_utils.calculate_norm());
    Aij_x.mult(m_utils.calculate_norm());
    scalar_fields_x.mult(m_utils.calculate_norm());

    // Convert to BSSN variables using the BSSN-CPT dictionary
    Aij_x.mult(-0.5);

    const auto &state_arrs          = state.arrays();
    const auto &hij_x_arrs          = hij_x.const_arrays();
    const auto &Aij_x_arrs          = Aij_x.const_arrays();
    const auto &scalar_field_x_arrs = scalar_fields_x.const_arrays();

    // Local copy so the kernel captures config by value, not via the host
    // `this` pointer
    const InflatonParameters d_params = params();

    amrex::ParallelFor(
        state,
        [=] AMREX_GPU_DEVICE(int bx, int i, int j, int k) noexcept
        {
            const amrex::IntVect iv_ds{i, j, k}; // coarse (state) index
            const amrex::IntVect iv{i * dn_ratio, j * dn_ratio,
                                    k * dn_ratio}; // fine (field) index

            if (iv_ds.min() >= 0 && iv_ds.max() < d_params.N)
            {
                const auto state_cell = state_arrs[bx];

                // Add scalar perturbations to the existing background values
                if (d_params.scalar_init)
                {
                    const auto scalar_x = scalar_field_x_arrs[bx];
                    state_cell(iv_ds, c_phi) +=
                        scalar_x(iv, static_cast<int>(BSSNFields::Phi));
                    state_cell(iv_ds, c_Pi) +=
                        scalar_x(iv, static_cast<int>(BSSNFields::Pi));
                    state_cell(iv_ds, c_chi) +=
                        scalar_x(iv, static_cast<int>(BSSNFields::Chi));
                    state_cell(iv_ds, c_K) +=
                        scalar_x(iv, static_cast<int>(BSSNFields::K));
                }

                // Add tensor perturbations to the existing background values
                if (d_params.tensor_init)
                {
                    const auto h_x = hij_x_arrs[bx];
                    const auto A_x = Aij_x_arrs[bx];
                    for (int comp = 0; comp < 6; comp++)
                    {
                        state_cell(iv_ds, c_h11 + comp) += h_x(iv, comp);
                        state_cell(iv_ds, c_A11 + comp) += A_x(iv, comp);
                    }
                }
            }
        });
    amrex::Gpu::streamSynchronize();
}

// Main initialisation routine
inline void ScalarTensorInit::init(amrex::MultiFab &state)
{
    BL_PROFILE("ScalarTensorInit::init");

    const int N      = params().N;
    const int N_fine = params().N_fine;

    // Derive MultiFab ingredients from state (configuration space)
    amrex::BoxArray sba            = state.boxArray();
    amrex::DistributionMapping sdm = state.DistributionMap();

    // If coarse graining is requested, set up the coarse grid Ns
    int n_eff = N;
    int dN    = 1;
    if (N_fine != N)
    {
        n_eff = N_fine;
        dN    = N_fine / N;
    }

    // Set up the problem domain in Fourier space
    // And impose that MPI ranks only slice along the i index (for Nyquist
    // conditions)
    amrex::IntVect domain_low(0, 0, 0);
    amrex::BoxArray xba = (N_fine != 0 ? sba.refine(dN) : sba);

    amrex::IntVect k_domain_high(n_eff / 2, n_eff - 1, n_eff - 1);
    amrex::Box k_domain(domain_low, k_domain_high);
    constexpr amrex::Array<bool, AMREX_SPACEDIM> slicing{true, false, false};
    amrex::BoxArray kba =
        decompose(k_domain, amrex::ParallelContext::NProcsAll(), slicing);
    amrex::DistributionMapping kdm(kba);

    // Set up the MFs to store the in/out data sets
    amrex::cMultiFab hij_k(kba, kdm, 6, 0);
    amrex::cMultiFab Aij_k(kba, kdm, 6, 0);
    amrex::cMultiFab scalar_fields_k(kba, kdm, 4, 0);
    amrex::MultiFab hij_x(xba, sdm, 6, 0);
    amrex::MultiFab Aij_x(xba, sdm, 6, 0);
    amrex::MultiFab scalar_fields_x(xba, sdm, 4, 0);

    hij_k.setVal(0.0);
    Aij_k.setVal(0.0);
    scalar_fields_k.setVal(0.0);
    hij_x.setVal(0.0);
    Aij_x.setVal(0.0);
    scalar_fields_x.setVal(0.0);

    // Construct the Fourier transform
    amrex::IntVect x_domain_high(n_eff - 1, n_eff - 1, n_eff - 1);
    amrex::Box x_domain(domain_low, x_domain_high);

    amrex::FFT::R2C<amrex::Real> tensor_fft(
        x_domain, amrex::FFT::Info().setBatchSize(hij_k.nComp()));

    amrex::FFT::R2C<amrex::Real> scalar_fft(
        x_domain, amrex::FFT::Info().setBatchSize(scalar_fields_k.nComp()));

    // Generate stochastic initial data
    generate_fourier_realisation(hij_k, Aij_k, scalar_fields_k);

    // Do the Fourier transform
    tensor_fft.backward(hij_k, hij_x);
    tensor_fft.backward(Aij_k, Aij_x);
    scalar_fft.backward(scalar_fields_k, scalar_fields_x);

    // Add perturbations to state
    add_perturbations_to_state(state, hij_x, Aij_x, scalar_fields_x, dN);
}

#endif /* SCALARTENSORINIT_IMPL_HPP_ */
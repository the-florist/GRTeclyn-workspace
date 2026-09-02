/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#if !defined(RANDOMFIELDINIT_HPP_)
#error "This file should only be included via RandomFieldInit.hpp"
#endif

#ifndef RANDOMFIELDINIT_IMPL_HPP_
#define RANDOMFIELDINIT_IMPL_HPP_

// Returns analytic power spectrum in modulus/argument form
AMREX_GPU_HOST_DEVICE inline amrex::GpuComplex<amrex::Real>
RandomFieldInit::calculate_mode_function(const InflationConfig &cfg,
                                         const amrex::Real km,
                                         const FieldType field_type,
                                         const WhichField which_field)
{
    // Deals with k=0 case, which is undefined if m=0
    if (km < InflationUtils::tolerance)
    {
        return amrex::GpuComplex<amrex::Real>{0., 0.};
    }

    // Stores modulus and argument
    amrex::Real ms_mag = 0.;
    amrex::Real ms_arg = 0.;

    amrex::Real kpr = km / cfg.H0;
    if (field_selector == WhichField::Amplitude) // Position mode funcion
    {
        ms_mag = sqrt((1.0 / km + cfg.H0 * cfg.H0 / pow(km, 3.)) / 2. /
                      pow(cfg.Mp, 2.));
        ms_arg =
            atan2((cos(kpr) + kpr * sin(kpr)), (kpr * cos(kpr) - sin(kpr)));
    }
    else // Velocity mode funcion
    {
        ms_mag = sqrt(km / 2. / pow(cfg.Mp, 2.));
        ms_arg = -atan2(cos(kpr), sin(kpr));
    }

    // Construct the mode function and return it
    amrex::GpuComplex<amrex::Real> ps(ms_mag * cos(ms_arg),
                                      ms_mag * sin(ms_arg));

    if (field_type == FieldType::Scalar)
    {
        ps /= std::sqrt(2. * cfg.epsilon_1); 
    }

    return ps;
}

// Turns analytic PS into GRF and applies window function if requested
AMREX_GPU_HOST_DEVICE inline amrex::GpuComplex<amrex::Real>
RandomFieldInit::calculate_random_field(const InflationConfig &cfg,
                                        const amrex::IntVect iv,
                                        const amrex::Real rand_amp,
                                        const amrex::Real rand_phase,
                                        const FieldType field_type,
                                        const WhichField which_field)
{
    amrex::GpuComplex<amrex::Real> value(0., 0.);
    amrex::Real kmag = cfg.get_kmag(iv);
    value = calculate_mode_function(cfg, kmag, field_type, which_field);

    // Add stochastic perturbations
    if (cfg.use_rand == 1)
    {
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
    }

    // Apply a window function if requested
    if (cfg.use_window == 1)
    {
        value *= cfg.window_function(kmag);
    }

    return value;
}

inline void
RandomFieldInit::generate_fourier_realisation(amrex::cMultiFab &hij_k,
                                              amrex::cMultiFab &Aij_k,
                                              amrex::cMultiFab &scalar_fields_k)
{
    // Test polarisation tensor orthonormality conditions
    if (inflt_methods.tensor_init && inflt_methods.test_normalisation)
    {
        inflt_methods.test_polarisation_normalisation(hij_k);
    }

    // Declare array to hold R and dR fields
    amrex::MultiFab R_dR_k(scalar_fields_k.boxArray(), 
                           scalar_fields_k.DistributionMap(),
                           2, 
                           scalar_fields_k.nGrowVect(), 
                           amrex::MFInfo(),
                           scalar_fields_k.Factory());
    R_dR_k.setVal(0.0);

    // Extract arrays before ParallelFor
    const auto &hij_k_arrs        = hij_k.arrays();
    const auto &Aij_k_arrs        = Aij_k.arrays();
    const auto &scalar_field_arrs = scalar_fields_k.arrays();
    const auto &R_dR_k_arrs       = R_dR_k.arrays();

    amrex::Print() << "Starting initial condition generation/read in...\n";

    // Slice to the POD base so the kernel captures config by value
    const InflationConfig cfg = inflt_methods;

    amrex::ParallelFor(
        hij_k,
        [=] AMREX_GPU_DEVICE(int bx, int i, int j, int k)
        {
            amrex::IntVect iv = {i, j, k};
            amrex::GpuArray<amrex::GpuComplex<amrex::Real>, 2>
                hs; // Amp mode fn
            amrex::GpuArray<amrex::GpuComplex<amrex::Real>, 2>
                As; // Vel mode fn

            // Uniform random draw, MPI/OpenMP safe
            const uint64_t cell_key =
                uint64_t(cfg.random_seed) *
                (645950ULL * uint64_t(iv[0]) +
                 520666ULL * uint64_t(cfg.invert_index_with_sign(iv[1])) +
                 767051ULL * uint64_t(cfg.invert_index_with_sign(iv[2])));

            // Initialise scalar sector (one random draw)
            if (cfg.scalar_init)
            {
                amrex::Real draw1 = InflationUtils::to_unit_open(
                    InflationUtils::splitmix64(cell_key + 0ULL));
                amrex::Real draw2 = InflationUtils::to_unit_open(
                    InflationUtils::splitmix64(cell_key + 1ULL));

                R_dR_k_arrs[bx](i, j, k, WhichField::Amplitude) = 
                    calculate_random_field(cfg, iv, draw1, draw2, 
                                           FieldType::Scalar,
                                           WhichField::Amplitude);
                
                R_dR_k_arrs[bx](i, j, k, WhichField::Velocity) = 
                    calculate_random_field(cfg, iv, draw1, draw2, 
                                           FieldType::Scalar,
                                           WhichField::Velocity);
            }

            // Initialise tensor sector (two random draws)
            if (cfg.tensor_init)
            {
                // Find the mode function realisation
                for (int p = 0; p < 2; p++)
                {
                    // One draw per polarisation field
                    amrex::Real draw1 =
                        InflationUtils::to_unit_open(InflationUtils::splitmix64(
                            cell_key + uint64_t(2 + 2 * p)));
                    amrex::Real draw2 =
                        InflationUtils::to_unit_open(InflationUtils::splitmix64(
                            cell_key + uint64_t(3 + 2 * p)));

                    hs[p] = calculate_random_field(cfg, iv, draw1, draw2,
                                                   FieldType::Tensor,
                                                   TensorField::Amplitude);

                    As[p] = calculate_random_field(cfg, iv, draw1, draw2,
                                                   FieldType::Tensor,
                                                   TensorField::Velocity);
                }

                // Construct polarisation tensors from basis vectors
                const auto [eplus, ecross] =
                    cfg.calculate_polarisation_tensors(iv);

                // Construct Fourier space tensors
                for (int l = 0; l < 3; l++)
                    for (int p = 0; p < 3; p++)
                    {
                        hij_k_arrs[bx](i, j, k, InflationUtils::lut[l][p]) =
                            (hs[0] * eplus(l, p) + hs[1] * ecross(l, p));
                        Aij_k_arrs[bx](i, j, k, InflationUtils::lut[l][p]) =
                            (As[0] * eplus(l, p) + As[1] * ecross(l, p));
                    }
            }
        });

    amrex::Gpu::streamSynchronize();

    // The R, dR -> BSSN functionals can go here, 
    // and should put the BSSN fields into the scalar_fields_k MF.
    // The scalar_fields_k MF is indexed by 0-3, 
    // and can be accessed using the BSSNFields enum, 
    // like in the assignment to state_cell below.
    // the slow-roll parameters can be accessed with cfg.epsilon_1 and 
    // cfg.epsilon_2, the Hubble parameter is cfg.H0 and Mp is cfg.Mp.
    // All possible parameters are found in the Config.hpp file.


    

    // Apply the DC and Nyquist symmetry conditions
    inflt_methods.apply_nyquist_conditions(hij_k);
    inflt_methods.apply_nyquist_conditions(Aij_k);
    inflt_methods.apply_nyquist_conditions(scalar_fields_k);
}

inline void RandomFieldInit::add_perturbations_to_state(
    amrex::MultiFab &state, amrex::MultiFab &hij_x, amrex::MultiFab &Aij_x,
    amrex::MultiFab &scalar_fields_x, const int dN)
{
    static_assert(c_h33 == c_h11 + 5 && c_A33 == c_A11 + 5,
                  "add_perturbations_to_state assumes c_h11..c_h33 and "
                  "c_A11..c_A33 are consecutive metric/A_ij components");

    // Apply normalisation into physical units
    hij_x.mult(inflt_methods.norm());
    Aij_x.mult(inflt_methods.norm());
    scalar_fields_x.mult(inflt_methods.norm());

    // Test that the resuling tensor perturbation field is trace-free
    TensorTests::Test_is_trace_free(hij_x);

    // Convert to BSSN variables using the BSSN-CPT dictionary
    Aij_x.mult(-0.5);

    const auto &state_arrs          = state.arrays();
    const auto &hij_x_arrs          = hij_x.const_arrays();
    const auto &Aij_x_arrs          = Aij_x.const_arrays();
    const auto &scalar_field_x_arrs = scalar_fields_x.const_arrays();

    // Slice to the POD base so the kernel captures config by value
    const InflationConfig cfg = inflt_methods;

    amrex::ParallelFor(
        state,
        [=] AMREX_GPU_DEVICE(int bx, int i, int j, int k) noexcept
        {
            const amrex::IntVect iv_ds{i, j, k}; // coarse (state) index
            const amrex::IntVect iv{i * dN, j * dN,
                                    k * dN}; // fine (field) index

            if (iv_ds.min() >= 0 && iv_ds.max() < cfg.N)
            {
                const auto state_cell = state_arrs[bx];

                // Add scalar perturbations to the existing background values
                if (cfg.scalar_init)
                {
                    const auto scalar_x       = scalar_field_x_arrs[bx];
                    state_cell(iv_ds, c_phi) += scalar_x(iv, BSSNFields::Phi);
                    state_cell(iv_ds, c_Pi)  += scalar_x(iv, BSSNFields::Pi);
                    state_cell(iv_ds, c_chi) += scalar_x(iv, BSSNFields::Chi);
                    state_cell(iv_ds, c_K)   += scalar_x(iv, BSSNFields::K);
                }

                // Add tensor perturbations to the existing background values
                if (cfg.tensor_init)
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
inline void RandomFieldInit::init(amrex::MultiFab &state)
{
    BL_PROFILE("RandomFieldInit::init");

    const int N      = inflt_methods.N;
    const int N_fine = inflt_methods.N_fine;

    // Derive MultiFab ingredients from state (configuration space)
    amrex::BoxArray sba            = state.boxArray();
    amrex::DistributionMapping sdm = state.DistributionMap();

    // If coarse graining is requested, set up the coarse grid Ns
    int Ni = N;
    int dN = 1;
    if (N_fine != N)
    {
        Ni = N_fine;
        dN = N_fine / N;
    }

    // Set up the problem domain in Fourier space
    // And impose that MPI ranks only slice along the i index (for Nyquist
    // conditions)
    amrex::IntVect domain_low(0, 0, 0);
    amrex::BoxArray xba = (N_fine != 0 ? sba.refine(dN) : sba);

    amrex::IntVect k_domain_high(Ni / 2, Ni - 1, Ni - 1);
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
    amrex::IntVect x_domain_high(Ni - 1, Ni - 1, Ni - 1);
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

#endif /* RANDOMFIELDINIT_IMPL_HPP_ */
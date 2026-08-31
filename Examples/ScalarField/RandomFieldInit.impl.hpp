/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */


#if !defined(RANDOMFIELDINIT_HPP_)
#error "This file should only be included via RandomFieldInit.hpp"
#endif

#ifndef RANDOMFIELDINIT_IMPL_HPP_
#define RANDOMFIELDINIT_IMPL_HPP_

// Estimate the loss of precision from adding a perturbation (the minimum
// absolute value in component comp of field) onto a background value bkgd.
// Errors if the perturbation is of the same order as, or larger than, bkgd.
inline amrex::Real RandomFieldInit::find_precision_loss(amrex::MultiFab &field, const int comp,
                                                        const amrex::Real bkgd)
{
    // 1. Initialize the reduction operator for a 'Minimum' operation
    amrex::ReduceOps<amrex::ReduceOpMin> reduce_op;
    amrex::ReduceData<amrex::Real> reduce_data(reduce_op);
    using ReduceTuple = typename decltype(reduce_data)::Type;

    // 2. Loop over the MultiFab (GPU and CPU safe)
    for (amrex::MFIter mfi(field); mfi.isValid(); ++mfi)
    {
        const amrex::Box &bx = mfi.fabbox();
        auto const &arr = field.array(mfi);

        reduce_op.eval(bx, reduce_data,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
            {
                // Return the absolute value of the cell
                return { std::abs(arr(i, j, k, comp)) };
            });
    }

    // 3. Extract the local minimum on this MPI rank
    amrex::Real min_abs_val = amrex::get<0>(reduce_data.value());

    // 4. Collective MPI reduction to find the global minimum across ranks
    amrex::ParallelDescriptor::ReduceRealMin(min_abs_val);

    int p_field = std::round(std::log10(min_abs_val));
    int p_bkgd = std::round(std::log10(std::abs(bkgd)));

    if (p_bkgd + p_field > 0)
    {
        amrex::Print() << bkgd << ", " << min_abs_val << "\n";
        amrex::Print() << p_bkgd << ", " << p_field << "\n";
        amrex::Error("RandomFieldInit::find_precision_loss, field may be non-perturbative.");
    }

    return pow(10., p_bkgd + p_field);
}

inline amrex::GpuComplex<amrex::Real> RandomFieldInit::find_in_stoiic(const amrex::Real km, const int field_indx,
                                                    const FieldType field_type)
{
    // Assume no average
    if (km == 0) { return amrex::GpuComplex<amrex::Real>{0., 0.}; }

    // Find the index where this k appears
    int spec_index = -1;
    for(int idx = 0; idx < m_params.init_k.size(); idx++)
    {
        if (std::abs(km - m_params.init_k[idx]) < InflationUtils::tolerance)
        {
            spec_index = idx;
            break;
        }
    }

    if (spec_index == -1)
    {
        amrex::Print() << km << "\n";
        amrex::Error("RandomFieldInit::find_in_stoiic, "
                     "the above k was not found in the STOIIC file.");
    }

    // Return the field at this k
    if (field_type == FieldType::Tensor)
    {
        return (amrex::GpuComplex<amrex::Real>{m_params.tensor_ps[2*field_indx][spec_index],
                                 m_params.tensor_ps[2*field_indx+1][spec_index]});
    }
    else
    {
        return (amrex::GpuComplex<amrex::Real>{m_params.scalar_ps[2*field_indx][spec_index],
                                 m_params.scalar_ps[2*field_indx+1][spec_index]});
    }
}

// Returns analytic power spectrum in modulus/argument form
inline amrex::GpuComplex<amrex::Real> RandomFieldInit::calculate_mode_function(const amrex::Real km, 
                                                             const int spec_indx)
{
    // Deals with k=0 case, which is undefined if m=0
    if (km < InflationUtils::tolerance) { return amrex::GpuComplex<amrex::Real>{0., 0.}; }
    
    // Stores modulus and argument 
    amrex::Real ms_mag = 0.;
    amrex::Real ms_arg = 0.;

    amrex::Real kpr = km/H0;
    if (spec_indx == 0) // Position mode funcion
    {
        ms_mag = sqrt((1.0/km + H0*H0/pow(km, 3.))/2./pow(m_params.Mp, 2.));
        ms_arg = atan2((cos(kpr) + kpr*sin(kpr)), (kpr*cos(kpr) - sin(kpr)));
    }
    else if (spec_indx == 1) // Velocity mode funcion
    {
        ms_mag = sqrt(km/2./pow(m_params.Mp, 2.));
        ms_arg = -atan2(cos(kpr), sin(kpr));
    }
    else { amrex::Error("RandomFieldInit::calculate_mode_function Value of spec_type not allowed."); }

    // Construct the mode function and return it
    amrex::GpuComplex<amrex::Real> ps(ms_mag * cos(ms_arg), ms_mag * sin(ms_arg));
    return ps;
}

// Turns analytic PS into GRF and applies window function if requested
inline amrex::GpuComplex<amrex::Real> RandomFieldInit::calculate_random_field(const amrex::IntVect iv, const int field_index,
                                                            const amrex::Real rand_amp, const amrex::Real rand_phase,
                                                            const FieldType field_type)
{
    amrex::GpuComplex<amrex::Real> value(0., 0.);

    // Find kmag with FFTW-style inversion on the last two indices
    amrex::Real kmag = m_params.get_kmag(iv);

    // Find the analytic power spectrum
    if (m_params.read_from_stoiic) 
    { 
        value = find_in_stoiic(kmag, field_index, field_type); 
    }
    else { value = calculate_mode_function(kmag, field_index); }

    // Add stochastic perturbations
    if (m_params.use_rand == 1)
    {
        BL_PROFILE("RandomFieldInit::calculate_random_field Random initialisation is used");

        // Make one random draw for the amplitude and phase individually
        amrex::Real rand_mod = sqrt(-2. * log(rand_amp)); // Rayleigh distribution about |h|
        amrex::Real rand_arg = 2. * M_PI * rand_phase;      // Uniform random phase

        // Multiply amplitude by Rayleigh draw
        value *= rand_mod;

        // Apply the random phase (assuming MS phase is accounted for)
        amrex::Real new_real = value.real() * cos(rand_arg) - value.imag() * sin(rand_arg);
        amrex::Real new_imag = value.real() * sin(rand_arg) + value.imag() * cos(rand_arg);
        amrex::GpuComplex<amrex::Real> new_value(new_real, new_imag);
	
        value = new_value;
    }

    // Apply a window function if requested
    if (m_params.use_window == 1)
    {
        BL_PROFILE("RandomFieldInit::calculate_random_field Window function is used")
        value *= m_params.window_function(kmag);
    }

    return value;
}

// Main initialisation routine
inline void RandomFieldInit::init(amrex::MultiFab &state)
{
    BL_PROFILE("RandomFieldInit::init");

    // Derive MultiFab ingredients from state (configuration space)
    amrex::BoxArray sba = state.boxArray();
    amrex::DistributionMapping sdm = state.DistributionMap();

    // If coarse graining is requested, set up the coarse grid Ns
    int Ni = m_params.N;
    int dN = 1;
    if(m_params.N_fine != m_params.N) 
    { 
        Ni = m_params.N_fine; 
        dN = m_params.N_fine / m_params.N; 
    }

    // Set up the problem domain in Fourier space
    // And impose that MPI ranks only slice along the i index (for Nyquist conditions)
    amrex::IntVect domain_low(0, 0, 0);
    amrex::BoxArray xba = (m_params.N_fine != 0 ? sba.refine(dN) : sba);

    amrex::IntVect k_domain_high(Ni/2, Ni-1, Ni-1);
    amrex::Box k_domain(domain_low, k_domain_high);
    constexpr amrex::Array<bool, AMREX_SPACEDIM> slicing{true, false, false};
    amrex::BoxArray kba = decompose(k_domain, amrex::ParallelContext::NProcsAll(), slicing);
    amrex::DistributionMapping kdm(kba);

    // Set up the MFs to store the in/out data sets
    amrex::cMultiFab hs_k(kba, kdm, 2, 0);
    amrex::cMultiFab As_k(kba, kdm, 2, 0);
    amrex::cMultiFab hij_k(kba, kdm, 6, 0);
    amrex::cMultiFab Aij_k(kba, kdm, 6, 0);
    amrex::MultiFab hij_x(xba, sdm, 6, 0);
    amrex::MultiFab Aij_x(xba, sdm, 6, 0);

    amrex::cMultiFab scalar_fields_k(kba, kdm, 4, 0);
    amrex::MultiFab scalar_fields_x(xba, sdm, 4, 0);

    hs_k.setVal(0.0);
    As_k.setVal(0.0);
    hij_k.setVal(0.0);
    Aij_k.setVal(0.0);
    hij_x.setVal(0.0);
    Aij_x.setVal(0.0);
    scalar_fields_k.setVal(0.0);
    scalar_fields_x.setVal(0.0);

    // One-time validation that the polarisation basis is orthonormal across
    // the Fourier grid (host-side; runs once here rather than per field point)
    if (m_params.tensor_init && m_params.test_normalisation)
    {
        m_params.test_polarisation_normalisation(hs_k);
    }

    // Construct the Fourier transform
    amrex::IntVect x_domain_high(Ni-1, Ni-1, Ni-1);
    amrex::Box x_domain(domain_low, x_domain_high);
    amrex::FFT::R2C<amrex::Real> tensor_fft(x_domain, amrex::FFT::Info().setBatchSize(hij_k.nComp()));
    amrex::FFT::R2C<amrex::Real> scalar_fft(x_domain, amrex::FFT::Info().setBatchSize(scalar_fields_k.nComp()));

    const auto &hs_arrs = hs_k.arrays();
    const auto &hij_k_arrs = hij_k.arrays();
    const auto &As_arrs = As_k.arrays();
    const auto &Aij_k_arrs = Aij_k.arrays();
    const auto &scalar_field_arrs = scalar_fields_k.arrays();

    amrex::Print() << "Starting initial condition generation/read in...\n";
    
    amrex::ParallelFor(hs_k,
        [=, this] AMREX_GPU_DEVICE (int bx, int i, int j, int k)
    {
        amrex::IntVect iv = {i, j, k};

        const uint64_t cell_key = uint64_t(m_params.random_seed)
            * (645950ULL * uint64_t(iv[0])
             + 520666ULL * uint64_t(m_params.invert_index_with_sign(iv[1]))
             + 767051ULL * uint64_t(m_params.invert_index_with_sign(iv[2]))
              );

        if (m_params.scalar_init)
        {
            amrex::Real draw1 = InflationUtils::to_unit_open(InflationUtils::splitmix64(cell_key + 0ULL));
            amrex::Real draw2 = InflationUtils::to_unit_open(InflationUtils::splitmix64(cell_key + 1ULL));
            for(int f=0; f<4; f++)
            {
                scalar_field_arrs[bx](i, j, k, f) = calculate_random_field(iv, f, draw1, draw2, FieldType::Scalar);
            }
        }

        if (m_params.tensor_init)
        {
            // Find the mode function realisation
            for(int p=0; p<2; p++)
            {
                amrex::Real draw1 = InflationUtils::to_unit_open(InflationUtils::splitmix64(cell_key + uint64_t(2 + 2*p)));
                amrex::Real draw2 = InflationUtils::to_unit_open(InflationUtils::splitmix64(cell_key + uint64_t(3 + 2*p)));

                hs_arrs[bx](i, j, k, p) = calculate_random_field(iv, 0, draw1, draw2, FieldType::Tensor);
                As_arrs[bx](i, j, k, p) = calculate_random_field(iv, 1, draw1, draw2, FieldType::Tensor);
            }
            
            // Construct polarisation tensors from basis vectors
            Tensor<2, amrex::Real> eplus = m_params.calculate_polarisation_tensor(iv, 0);
            Tensor<2, amrex::Real> ecross = m_params.calculate_polarisation_tensor(iv, 1);

            // Find basis tensors and initial tensor realisation
            for (int l=0; l<3; l++) for (int p=0; p<3; p++)
            {
                hij_k_arrs[bx](i, j, k, InflationUtils::lut[l][p]) = (hs_arrs[bx](i, j, k, 0) * eplus[l][p] 
                                                                    + hs_arrs[bx](i, j, k, 1) * ecross[l][p]);
                Aij_k_arrs[bx](i, j, k, InflationUtils::lut[l][p]) = (As_arrs[bx](i, j, k, 0) * eplus[l][p] 
                                                                    + As_arrs[bx](i, j, k, 1) * ecross[l][p]);
            }
        }
    });

    amrex::Gpu::streamSynchronize();

    // Apply the DC and Nyquist symmetry conditions
    m_params.apply_nyquist_conditions(hij_k);
    m_params.apply_nyquist_conditions(Aij_k);
    m_params.apply_nyquist_conditions(scalar_fields_k);

    // Do the Fourier transform
    tensor_fft.backward(hij_k, hij_x);
    tensor_fft.backward(Aij_k, Aij_x);
    scalar_fft.backward(scalar_fields_k, scalar_fields_x);

    // Apply normalisation into physical units
    hij_x.mult(m_params.norm());
    Aij_x.mult(m_params.norm());
    scalar_fields_x.mult(m_params.norm());

    // Check the scalar perturbations remain perturbative on the background
    if (m_params.scalar_init)
    {
        amrex::Print() << "RandomFieldInit::init, Precision lost in phi is ";
        amrex::Print() << find_precision_loss(scalar_fields_x, 0, phi0) << "\n";
        amrex::Print() << "RandomFieldInit::init, Precision lost in chi is ";
        amrex::Print() << find_precision_loss(scalar_fields_x, 2, 1.0) << "\n";
    }

    // Test that the resuling tensor perturbation field is trace-free
    TensorTests::Test_is_trace_free(hij_x);

    // Convert to BSSN variables using the BSSN-CPT dictionary
    Aij_x.mult(-0.5);

    const auto &state_arrs = state.arrays();
    const auto &hij_x_arrs = hij_x.arrays();
    const auto &Aij_x_arrs = Aij_x.arrays();
    const auto &scalar_field_x_arrs = scalar_fields_x.arrays();

    amrex::ParallelFor(hij_x,
        [=, this] AMREX_GPU_DEVICE (int bx, int i, int j, int k) noexcept
    {
        const amrex::IntVect iv_ds{i, j, k};
        const amrex::IntVect iv{i * dN, j * dN, k * dN};

        if (iv_ds.min() >= 0 && iv_ds.max() < m_params.N)
        {
            // Add scalar perturbations to the existing background values
            if (m_params.scalar_init)
            {
                state_arrs[bx](iv_ds, c_phi) += scalar_field_x_arrs[bx](i, j, k, 0);
                state_arrs[bx](iv_ds, c_Pi) += scalar_field_x_arrs[bx](i, j, k, 1);
                state_arrs[bx](iv_ds, c_chi) += scalar_field_x_arrs[bx](i, j, k, 2);
                state_arrs[bx](iv_ds, c_K) += scalar_field_x_arrs[bx](i, j, k, 3);
            }

            // Set the entire tensor object here
            if (m_params.tensor_init)
            {
                state_arrs[bx](iv_ds, c_h11) += hij_x_arrs[bx](i, j, k, InflationUtils::lut[0][0]);
                state_arrs[bx](iv_ds, c_h12) += hij_x_arrs[bx](i, j, k, InflationUtils::lut[0][1]);
                state_arrs[bx](iv_ds, c_h13) += hij_x_arrs[bx](i, j, k, InflationUtils::lut[0][2]);
                state_arrs[bx](iv_ds, c_h22) += hij_x_arrs[bx](i, j, k, InflationUtils::lut[1][1]);
                state_arrs[bx](iv_ds, c_h23) += hij_x_arrs[bx](i, j, k, InflationUtils::lut[1][2]);
                state_arrs[bx](iv_ds, c_h33) += hij_x_arrs[bx](i, j, k, InflationUtils::lut[2][2]);

                state_arrs[bx](iv_ds, c_A11) += Aij_x_arrs[bx](i, j, k, InflationUtils::lut[0][0]);
                state_arrs[bx](iv_ds, c_A12) += Aij_x_arrs[bx](i, j, k, InflationUtils::lut[0][1]);
                state_arrs[bx](iv_ds, c_A13) += Aij_x_arrs[bx](i, j, k, InflationUtils::lut[0][2]);
                state_arrs[bx](iv_ds, c_A22) += Aij_x_arrs[bx](i, j, k, InflationUtils::lut[1][1]);
                state_arrs[bx](iv_ds, c_A23) += Aij_x_arrs[bx](i, j, k, InflationUtils::lut[1][2]);
                state_arrs[bx](iv_ds, c_A33) += Aij_x_arrs[bx](i, j, k, InflationUtils::lut[2][2]);
            }
        }
    });
    amrex::Gpu::streamSynchronize();
}

#endif /* RANDOMFIELDINIT_IMPL_HPP_ */
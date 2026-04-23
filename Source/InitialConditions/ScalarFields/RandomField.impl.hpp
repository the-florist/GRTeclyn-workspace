/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */


#if !defined(RANDOMFIELD_HPP_)
#error "This file should only be included via RandomField.hpp"
#endif

#ifndef RANDOMFIELD_IMPL_HPP_
#define RANDOMFIELD_IMPL_HPP_

/****
    Small functions (less than 10 lines)
****/

// Nyquist condition
inline int RandomField::flip_index(const int indx) { return std::abs(N - indx); }

// Nyquist condition and calculation of kmag
inline int RandomField::invert_index(const int indx) { return (int)(N/2 - std::abs(N/2 - indx)); }

// For calculation of polarisation tensors
inline int RandomField::invert_index_with_sign(const int indx) 
{ 
    if(indx <= N/2) { return indx; }
    else { return std::abs(N/2 - indx) - N/2; }
}

// Find the magnitude of the Fourier wavevector at this point
inline Real RandomField::get_kmag(IntVect iv)
{
    const int i = iv[0];
    const int j = invert_index(iv[1]);
    const int k = invert_index(iv[2]);
    return std::sqrt(i*i + j*j + k*k) * 2. * M_PI / m_params.L;
}

// Ensures no calculation on ghost cells
inline bool RandomField::is_ghost_index(const IntVect vector)
{
    bool ret = false;
    for(int d=0; d<3; d++) 
    { 
        if(vector[d] < 0 || vector[d] > N-1) { ret = true; }
    }
    return ret;
}

// Makes subdirectories in data/
inline std::string RandomField::make_subdirectory(const std::string base, const std::string dir, const int is_first_step)
{
    std::string new_path = base+"../"+dir+"/";
    if(is_first_step)
    {
        if (FilesystemTools::directory_exists(base)) { FilesystemTools::mkdir_recursive(new_path); }
        else 
        { 
            Print() << "RandomField::make_subdirectory, Directory creation failed for " << new_path << "\n";
            Error("RandomField::extract Data directory has not been created."); 
        }
    }
    return new_path;
}

// Creates a custom data file layout 
inline void RandomField::assign_statistics_data(Vector<std::string> &header_storage, const std::string name, 
                            Vector<Real> &data_storage, const Vector<Real> data, const int component, const int num_comps,
                            const Vector<int>::const_iterator itr, const Vector<int>::const_iterator start, 
                            const int is_first_step)
{
    int loc = component + num_comps*(itr - start);
    if(is_first_step) 
    { 
        header_storage[loc] =  name; 
    }
    data_storage[loc] = data[component];
}

// Written by Gemini, edited by me
inline Real RandomField::find_precision_loss(MultiFab &field, int comp, Real bkgd)
{
    // 1. Initialize the reduction operator for a 'Minimum' operation
    amrex::ReduceOps<amrex::ReduceOpMin> reduce_op;
    amrex::ReduceData<amrex::Real> reduce_data(reduce_op);
    using ReduceTuple = typename decltype(reduce_data)::Type;

    // 2. Loop over the MultiFab (GPU and CPU safe)
    for (amrex::MFIter mfi(field); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.fabbox();
        auto const& arr = field.array(mfi);

        reduce_op.eval(bx, reduce_data,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple {
                // Return the absolute value of the cell
                return { std::abs(arr(i, j, k, comp)) };
            });
    }

    // 3. Extract the local minimum on this MPI rank
    amrex::Real min_abs_val = amrex::get<0>(reduce_data.value());

    // 4. Perform a collective MPI reduction to find the global minimum across all processors
    amrex::ParallelDescriptor::ReduceRealMin(min_abs_val);

    int p_field = std::round(std::log10(min_abs_val));
    int p_bkgd = std::round(std::log10(std::abs(bkgd)));

    if (p_bkgd + p_field > 0)
    {
        Print() << bkgd << ", " << min_abs_val << "\n";
        Print() << p_bkgd << ", " << p_field << "\n";
        Error("RandomField::find_precision_loss, field may be non-perturbative.");
    }

    return pow(10., p_bkgd + p_field);
}

/****
    Tests
****/

// Test that the input tensor field (config space) is trace free (global)
inline void RandomField::Test_is_trace_free(MultiFab &field)
{
    if (field.nComp() != 6)
    {
        Error("RandomField::Test_is_trace_free, input field is not a tensor field");
    }

    for (MFIter mfi(field); mfi.isValid(); ++mfi) 
    {
        Array4<Real> const& field_ptr = field.array(mfi);
        const Box& bx = mfi.fabbox();

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            IntVect iv{i, j, k};
            Real sum = 0.;

            for(int l=0; l<3; l++)
            {
                sum += field_ptr(i, j, k, lut[l][l]);
            }

            if(std::abs(sum) > tolerance)
            {
                Print() << iv << ": " << sum << "\n";
                Error("RandomField::Test_is_trace_free Trace-free test failed here.");
            }
        });
    }
    
}

// Test that the input vectors are orthonormal (local)
inline void RandomField::Test_vector_orthonorm(const IntVect iv, const Vector<Real> mhat, 
                                                                 const Vector<Real> nhat)
{
    // Confirm basis vectors are orthonormal
    if (iv != IntVect{0, 0, 0})
    {
        Real dot1 = 0.;
        Real dot2 = 0.;
        Real cross1 = 0.;
        for(int l=0; l<3; l++)
        {
            dot1 += mhat[l] * mhat[l];
            dot2 += nhat[l] * nhat[l];
            cross1 += mhat[l] * nhat[l];
        }

        if(std::abs(dot1 - 1.) > tolerance || std::abs(dot2 - 1.) > tolerance || cross1 > tolerance)
        {
            Print() << "Location: " << iv << "\n";
            Print() << "Dot products: " << dot1 << ", " << dot2 << "\n";
            Print() << "Cross product: " << cross1 << "\n";
            Print() << "Vectors: \n";
            for(int l=0; l<3; l++)
            {
                Print() << l << ", " << mhat[l] << ", " << nhat[l] << "\n";
            }
            Error("RandomField::Test_vector_orthonorm: Basis vectors are not orthonormal here");
        }
    }
}

// Test that the input basis tensors, and their rotated counterparts, are orthonormal
inline void RandomField::Test_polarisation_tensor_orthonorm(const IntVect iv, const Tensor<2, Real> eplus,
                                                                              const Tensor<2, Real> ecross)
{
    Vector<Real> conds(3, 0.);

    for (int l=0; l<3; l++) for (int p=0; p<3; p++)
    {
        conds[0] += eplus[l][p] * eplus[l][p];
        conds[1] += eplus[l][p] * ecross[l][p];
        conds[2] += ecross[l][p] * ecross[l][p];
    }

    if(iv != IntVect{0, 0, 0})
    {
        bool plc = (std::abs(conds[0] / 2. - 1.) > tolerance 
                    || std::abs(conds[2] / 2. - 1.) > tolerance);
        bool ppc = (std::abs(conds[1]) > tolerance);
        if (plc || ppc)
        {
            Print() << "---------\nLocation: " << iv << "\n";
            Print() << "Dot products: " << conds[0] / 2. << ", " << conds[2] / 2. << "\n";
            Print() << "Cross product: " << conds[1] << "\n";
            Print() << "Base tensor components: \n";
            for (int l=0; l<3; l++) for (int p=0; p<3; p++)
            {
                Print() << l << ", " << p << ": ";
                Print() << eplus[l][p] << ", " << ecross[l][p] << "\n";
            }
            Error("RandomField::Test_polarisation_tensor_orthonorm: polarisation tensors are not orthonormal here");
        }
    }
}

// Written by Gemini
inline Real RandomField::calculate_total_power(const cMultiFab& fk) 
{
    // 1. Set up the parallel reduction operation (Sum)
    ReduceOps<ReduceOpSum> reduce_op;
    ReduceData<Real> reduce_data(reduce_op);
    using ReduceTuple = typename decltype(reduce_data)::Type;

    // 2. Loop over the grids owned by this MPI rank
    for (MFIter mfi(fk); mfi.isValid(); ++mfi) 
    {
        const Box& bx = mfi.fabbox();
        auto const& arr = fk.const_array(mfi);

        // 3. Execute the parallel reduction on the CPU/GPU
        reduce_op.eval(bx, reduce_data,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
            {
                // Get the real and imaginary parts
                Real re = arr(i,j,k).real();
                Real im = arr(i,j,k).imag();

                Real pow = re * re + im * im;
                // Multiply by 2 in most of the bulk,
                // to account for Hermitian modes.
                if (i != 0 && i != N/2) { pow *= 2.; }
                
                return pow;
            });
    }

    // 4. Extract the aggregated sum for this specific MPI rank
    ReduceTuple hv = reduce_data.value();
    Real total_local_power = get<0>(hv);

    // 5. Perform an MPI reduction to sum across all compute nodes
    ParallelDescriptor::ReduceRealSum(total_local_power);

    return total_local_power;
}

inline void RandomField::Test_Parsevals_thm(const MultiFab &hx, const cMultiFab &hk)
{
    Real xsum = std::pow(hx.norm2(), 2.);
    xsum /= std::pow(N, 3.);

    Real ksum = calculate_total_power(hk);

    int p = std::round(std::log10((ksum + ksum) / 2.));
    Real tol = tolerance * std::pow(10., p);

    if (std::abs(xsum - ksum) > tol)
    {
        Print() << "Order of magnitude of sum: " << p << "\n";
        Print() << "Tolerance: " << tol << "\n";
        Print() << "Stdev (x): " << xsum << "\n";
        Print() << "Integrated power (k): " << ksum << "\n";
        Print() << "Ratio: " << ksum / xsum << "\n";
        Print() << "Difference: " << std::abs(ksum - xsum) << "\n";
        Error("RandomField::Test_Parsevals_thm, Parseval's theorem fails here.");
    }
}

/****
    Initialisation routines
****/

// Written by Gemini
inline Real RandomField::get_spatial_random(int i, int j, int k, int comp, int seed)
{
    // 1. Create a unique 64-bit identifier for this exact cell, component, and seed.
    // The large prime numbers help prevent spatial artifacts (striping) before mixing.
    uint64_t state = (uint64_t(i) * 73856093ULL) ^
                     (uint64_t(j) * 19349663ULL) ^
                     (uint64_t(k) * 83492791ULL) ^
                     (uint64_t(comp) * 23145671ULL) ^
                     (uint64_t(seed));

    // 2. High-quality bit mixer (based on the SplitMix64 algorithm)
    // This violently scrambles the bits so neighboring cells have no correlation.
    state ^= state >> 30;
    state *= 0xbf58476d1ce4e5b9ULL;
    state ^= state >> 27;
    state *= 0x94d049bb133111ebULL;
    state ^= state >> 31;

    // 3. Convert the scrambled 64-bit integer into a double precision float in [0.0, 1.0)
    // 0x1.0p-53 is a fast hex-float representation of 2^-53.
    return (state >> 11) * 0x1.0p-53;
}

// Returns analytic power spectrum in modulus/argument form
inline GpuComplex<Real> RandomField::calculate_mode_function(const Real km, const int spec_indx)
{
    // Deals with k=0 case, which is undefined if m=0
    if(km < 1.e-23) { return 0.; }
    
    // Stores modulus and argument 
    Real ms_mag = 0.;
    Real ms_arg = 0.;

    Real kpr = km/H0;
    if (spec_indx == 0) // Position mode funcion
    {
        ms_mag = sqrt((1.0/km + H0*H0/pow(km, 3.))/2.);
        ms_arg = atan2((cos(kpr) + kpr*sin(kpr)), (kpr*cos(kpr) - sin(kpr)));
    }
    else if (spec_indx == 1) // Velocity mode funcion
    {
        ms_mag = sqrt(km/2.);
        ms_arg = -atan2(cos(kpr), sin(kpr));
    }
    else { Error("RandomField::calculate_mode_function Value of spec_type not allowed."); }

    // Construct the mode function and return it
    GpuComplex<Real> ps(ms_mag * cos(ms_arg), ms_mag * sin(ms_arg));
    return ps;
}

inline GpuComplex<Real> RandomField::find_in_stoiic(const Real km, const int field_indx, std::string field_type)
{
    GpuComplex<Real> zero(0., 0.);
    if(km == 0) { return zero; }

    int spec_index;
    for(int idx = 0; idx < m_params.init_k.size(); idx++)
    {
        if(std::abs(km - m_params.init_k[idx]) < 1e-10) { spec_index = idx; break; }
        else if (idx == m_params.init_k.size() - 1) 
        { 
            Print() << km << "\n"; 
            Error("RandomField::find_in_stoiic, The above k was not found in the STOIIC file."); 
        }
    }

    if(field_type == "tensor")
    {
        return GpuComplex<Real>{m_params.tensor_ps[2*field_indx][spec_index], m_params.tensor_ps[2*field_indx+1][spec_index]};
    }
    else if(field_type == "scalar")
    {
        return GpuComplex<Real>{m_params.scalar_ps[2*field_indx][spec_index], m_params.scalar_ps[2*field_indx+1][spec_index]};
    }
    else { Error("RandomField::find_in_stoiic field cannot be found."); return GpuComplex<Real>{0., 0.}; }
}

// Turns analytic PS into GRF and applies window function if requested
inline GpuComplex<Real> RandomField::calculate_random_field(const IntVect iv, const int field_index, 
                                                            const Real rand_amp, const Real rand_phase, 
                                                            std::string field_type)
{
    GpuComplex<Real> value(0., 0.);

    // Find kmag with FFTW-style inversion on the last two indices
    Real kmag = get_kmag(iv);

    // Find the analytic power spectrum
    if(m_params.read_from_stoiic) { value = find_in_stoiic(kmag, field_index, field_type); }
    else { value = calculate_mode_function(kmag, field_index); }

    // Add stochastic perturbations
    if(m_params.use_rand == 1)
    {
        BL_PROFILE("RandomField::calculate_random_field Random initialisation is used");

        // Make one random draw for the amplitude and phase individually
        Real rand_mod = sqrt(-2. * log(rand_amp)); // Rayleigh distribution about |h|
        Real rand_arg = 2. * M_PI * rand_phase;      // Uniform random phase

        // Multiply amplitude by Rayleigh draw
        value *= rand_mod;

        // Apply the random phase (assuming MS phase is accounted for)
        Real new_real = value.real() * cos(rand_arg) - value.imag() * sin(rand_arg);
        Real new_imag = value.real() * sin(rand_arg) + value.imag() * cos(rand_arg);
        GpuComplex<Real> new_value(new_real, new_imag);
	
        value = new_value;
    }

    // Apply a window function if requested
    if(m_params.use_window == 1) 
    { 
        BL_PROFILE("RandomField::calculate_random_field Window function is used")
        Real ks = (m_params.N_coarse != 0 ? 
                   std::sqrt(3.) * m_params.N_coarse * M_PI / m_params.L / 5. / 2. :
                   std::sqrt(3.) * N * M_PI / m_params.L / 5. / 2.);
        //m_params.kstar * 2. * M_PI/m_params.L;
        Real Dt = m_params.L/m_params.Delta;
        value *= 0.5 * (1.0 - tanh(Dt * (kmag - ks))); 
    }

    return value;
}

// Calculates basis vectors required for polarisation tensors
inline Vector<Real> RandomField::calculate_basis_vector(const IntVect iv, const int which_vector)
{
    // Hermitian symmetry inversion on j and k, with sign on the last two indices.
    // (!!) The FT implemented in AMReX symmetrises across the i index.
    const Real i = static_cast<Real>(iv[0]);
    const Real j = static_cast<Real>(invert_index_with_sign(iv[1]));
    const Real k = static_cast<Real>(invert_index_with_sign(iv[2]));

    Vector<Real> mhat(3, 0.);
    Vector<Real> nhat(3, 0.);

    // Skip the 0 mode, as tensors have no average
    if (iv == IntVect{0, 0, 0}) { ; }

    else if (i > 0.) 
    {
        if (k == 0. && j == 0.) 
        { 
            mhat = Vector<Real>{0., 1., 0.};
            nhat = Vector<Real>{0., 0., 1.}; 
        }

        else 
        { 
            Real norm = sqrt((i*i + j*j) * (i*i + j*j + k*k));
            mhat = Vector<Real>{j/sqrt(i*i + j*j), -i/sqrt(i*i + j*j), 0.}; 
            nhat = Vector<Real>{(i*k) / norm, (j*k) / norm, -(i*i + j*j) / norm}; 
        }
    }

    else if (std::abs(j) > 0) 
    { 
        if(k == 0.)
        {
            mhat = Vector<Real>{0., 0., 1.};
            nhat = Vector<Real>{1., 0., 0.};
        }

        else
        {
            mhat = Vector<Real>{-1., 0., 0.};
            nhat = Vector<Real>{0., -k / sqrt(j*j + k*k), j / sqrt(j*j + k*k)};
        }
    }

    else if (std::abs(k) > 0) 
    { 
        mhat = Vector<Real>{1., 0., 0.};
        nhat = Vector<Real>{0., 1., 0.};
    }

    else 
    {
        Error("RandomField::calculate_polarisation_tensors Part of Fourier grid not covered.");
    }

    if (m_params.alpha != 0)
    {
        Real a = m_params.alpha * (M_PI) / 180.;
        Vector<Real> mp(3), np(3);
        for(int l=0; l<3; l++)
        {
            mp[l] = cos(a) * mhat[l] + sin(a) * nhat[l];
            np[l] = -sin(a) * mhat[l] + cos(a) * nhat[l];
        }

        if(which_vector == 0) { return mp; }
        else if(which_vector == 1) { return np; }
        else 
        { 
            Error("RandomField::calculate_basis_vector Incompatable vector type."); 
            return Vector<Real>{0,0,0}; 
        }
    }

    else
    {
        if(which_vector == 0) { return mhat; }
        else if(which_vector == 1) { return nhat; }
        else 
        { 
            Error("RandomField::calculate_basis_vector Incompatable vector type."); 
            return Vector<Real>{0,0,0}; 
        }
    }
}

// Applies above Nyquist conditions to a given MF
inline void RandomField::apply_nyquist_conditions(cMultiFab &field)
{
    int nc = field.nComp();
    for (MFIter mfi(field); mfi.isValid(); ++mfi) 
    {
        // The geometry for this MPI rank
        const Box& bx = mfi.fabbox();
        Array4<GpuComplex<Real>> const& field_ptr = field.array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            IntVect iv = {i, j, k};

            if ((i == 0 || i == N/2) && (j == 0 || j == N/2) && (k == 0 || k == N/2))
            {
                for(int comp = 0; comp < nc; comp++)
                {
                    GpuComplex<Real> temp(field_ptr(i, j, k, comp).real(), 0.);
                    field_ptr(i, j, k, comp) = temp;
                }
            }

            else if (i==0 || i==N/2) 
            {
                if((k > N/2 && j == N/2) || (k == 0 && j > N/2) ||
                    (k > N/2 && j == 0) || (k == N/2 && j > N/2))
                {
                    for(int comp = 0; comp < nc; comp++) 
                    {
                        GpuComplex<Real> temp(field_ptr(i, invert_index(j), invert_index(k), comp).real(), 
                                                -field_ptr(i, invert_index(j), invert_index(k), comp).imag());
                        field_ptr(i, j, k, comp) = temp;
                    }
                }
                
                else if(j > N/2)
                {
                    for(int comp = 0; comp < nc; comp++) 
                    {
                        GpuComplex<Real> temp(field_ptr(i, invert_index(j), flip_index(k), comp).real(), 
                                                -field_ptr(i, invert_index(j), flip_index(k), comp).imag());
                        field_ptr(i, j, k, comp) = temp;
                    }
                }
            }
        });
    }
}

// Main initialisation routine
inline void RandomField::init(amrex::MultiFab &state)
{
    BL_PROFILE("RandomField::init");

    // Derive MultiFab ingredients from state (configuration space)
    BoxArray sba = state.boxArray();
    DistributionMapping sdm = state.DistributionMap();

    // Set up the problem domain in Fourier space
    // And impose that MPI ranks only slice along the i index (for Nyquist conditions)
    IntVect domain_low(0, 0, 0);
    IntVect k_domain_high(N/2, N-1, N-1);
    Box k_domain(domain_low, k_domain_high);
    Array< bool, AMREX_SPACEDIM > const &slicing{true, false, false};
    BoxArray kba = decompose(k_domain, ParallelContext::NProcsAll(), slicing);
    DistributionMapping kdm(kba);

    // Set up the MFs to store the in/out data sets
    cMultiFab hs_k(kba, kdm, 2, 0);
    cMultiFab As_k(kba, kdm, 2, 0);
    cMultiFab hij_k(kba, kdm, 6, 0);
    cMultiFab Aij_k(kba, kdm, 6, 0);

    MultiFab hij_x(sba, sdm, 6, 0);
    MultiFab Aij_x(sba, sdm, 6, 0);

    cMultiFab scalar_fields_k(kba, kdm, 4, 0);
    MultiFab scalar_fields_x(sba, sdm, 4, 0);

    hs_k.setVal(0.0);
    As_k.setVal(0.0);
    hij_k.setVal(0.0);
    Aij_k.setVal(0.0);
    hij_x.setVal(0.0);
    Aij_x.setVal(0.0);
    scalar_fields_k.setVal(0.0);
    scalar_fields_x.setVal(0.0);

    // Construct the Fourier transform
    IntVect x_domain_high(N-1, N-1, N-1);
    Box x_domain(domain_low, x_domain_high);
    FFT::R2C<Real> tensor_fft(x_domain, FFT::Info().setBatchSize(hij_k.nComp()));
    FFT::R2C<Real> scalar_fft(x_domain, FFT::Info().setBatchSize(scalar_fields_k.nComp()));

    Print() << "RandomField::init, Starting initial condition generation/read in...\n";
    for (MFIter mfi(hs_k); mfi.isValid(); ++mfi) 
    {
        // Define the domain on this MPI rank
        const Box& bx = mfi.fabbox();

        // Make a pointer to the mode functions at this MPI box
        Array4<GpuComplex<Real>> const& hs_ptr = hs_k.array(mfi);
        Array4<GpuComplex<Real>> const& hij_ptr = hij_k.array(mfi);

        Array4<GpuComplex<Real>> const& As_ptr = As_k.array(mfi);
        Array4<GpuComplex<Real>> const& Aij_ptr = Aij_k.array(mfi);

        Array4<GpuComplex<Real>> const& scalar_fields_ptr = scalar_fields_k.array(mfi);

        // Loop to create mode functions, then hij(k) and Aij(k)
        amrex::ParallelForRNG(bx, 
            [=] AMREX_GPU_DEVICE (int i, int j, int k, 
                                  RandomEngine const& engine) noexcept
        {
            IntVect iv = {i, j, k};

            if(m_params.scalar_init)
            {
                Real draw1 = amrex::Random(engine);
                            /*get_spatial_random(i, invert_index_with_sign(j), 
                                                   invert_index_with_sign(k), 
                                                   4, m_params.random_seed);*/
                                                   
                Real draw2 = amrex::Random(engine);
                            /*get_spatial_random(i, invert_index_with_sign(j), 
                                                   invert_index_with_sign(k), 
                                                   5, m_params.random_seed);*/

                AllPrintToFile("dump/amrex-rand-check") << iv[0] << ", ";
                AllPrintToFile("dump/amrex-rand-check") << iv[1] << ", ";
                AllPrintToFile("dump/amrex-rand-check") << iv[2] << ", ";
                AllPrintToFile("dump/amrex-rand-check").SetPrecision(12) << draw1 << ", ";
                AllPrintToFile("dump/amrex-rand-check").SetPrecision(12) << draw2 << "\n";
                
                for(int f=0; f<4; f++)
                {
                    scalar_fields_ptr(i, j, k, f) = calculate_random_field(iv, f, draw1, draw2, "scalar");
                    if (m_params.A != 1.0)
                    {
                        scalar_fields_ptr(i, j, k, f) *= m_params.A;
                    }
                }
            }

            if(m_params.tensor_init)
            {
                // Find the mode function realisation
                for(int p=0; p<2; p++)
                {
                    Real draw1 = amrex::Random(engine);
                    //get_spatial_random(i, invert_index_with_sign(j), invert_index_with_sign(k), 2*p, m_params.random_seed);
                    Real draw2 = amrex::Random(engine);
                    //get_spatial_random(i, invert_index_with_sign(j), invert_index_with_sign(k), 2*p+1, m_params.random_seed);

                    hs_ptr(i, j, k, p) = calculate_random_field(iv, 0, draw1, draw2, "tensor");
                    As_ptr(i, j, k, p) = calculate_random_field(iv, 1, draw1, draw2, "tensor");
                }
                
                // Find basis vectors
                Vector<Real> mhat = calculate_basis_vector(iv, 0);;
                Vector<Real> nhat = calculate_basis_vector(iv, 1);
                Test_vector_orthonorm(iv, mhat, nhat);

                // Construct polarisation tensors from basis vectors
                Tensor<2, Real> eplus, ecross; 

                // Find basis tensors and initial tensor realisation
                for (int l=0; l<3; l++) for (int p=0; p<3; p++)
                {
                    // Assemble the polarisation tensors
                    eplus[l][p] = mhat[l]*mhat[p] - nhat[l]*nhat[p];
                    ecross[l][p] = mhat[l]*nhat[p] + nhat[l]*mhat[p];

                    hij_ptr(i, j, k, lut[l][p]) = hs_ptr(i, j, k, 0) * eplus[l][p] + hs_ptr(i, j, k, 1) * ecross[l][p];
                    Aij_ptr(i, j, k, lut[l][p]) = As_ptr(i, j, k, 0) * eplus[l][p] + As_ptr(i, j, k, 1) * ecross[l][p];
                }

                if (m_params.alpha != 0) { Test_polarisation_tensor_orthonorm(iv, eplus, ecross); }
            }
        });
    }

    // Apply the DC and Nyquist symmetry conditions
    apply_nyquist_conditions(hij_k);
    apply_nyquist_conditions(Aij_k);
    apply_nyquist_conditions(scalar_fields_k);

    // Do the Fourier transform
    tensor_fft.backward(hij_k, hij_x);
    tensor_fft.backward(Aij_k, Aij_x);
    scalar_fft.backward(scalar_fields_k, scalar_fields_x);

    // Apply normalisation into physical units
    hij_x.mult(norm * std::pow(N, -3./2.));
    Aij_x.mult(norm * std::pow(N, -3./2.));
    scalar_fields_x.mult(norm * std::pow(N, -3./2.));

    Print() << "RandomField::init, Precision lost in phi is ";
    Print() << find_precision_loss(scalar_fields_x, 0, phi0) << "\n";
    Print() << "RandomField::init, Precision lost in chi is ";
    Print() << find_precision_loss(scalar_fields_x, 2, 1.0) << "\n";

    // Test that the resuling tensor perturbation field is trace-free
    Test_is_trace_free(hij_x);

    // Convert to BSSN variables using the BSSN-CPT dictionary
    for (int l=0; l<3; l++) { hij_x.plus(1., lut[l][l], 1); }
    Aij_x.mult(-0.5);

    // Put these initial conditions into the state MF
    for (MFIter mfi(hij_x); mfi.isValid(); ++mfi) 
    {
        const Box& bx = mfi.fabbox();
        Array4<Real> const& state_ptr = state.array(mfi);
        Array4<Real> const& hij_ptr = hij_x.array(mfi);
        Array4<Real> const& Aij_ptr = Aij_x.array(mfi);
        Array4<Real> const& scalar_ptr = scalar_fields_x.array(mfi);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            const IntVect iv{i, j, k};
            bool in_ghost_index = is_ghost_index(iv);
            if(!in_ghost_index)
            {
                // Add scalar perturbations to the existing background values
                if(m_params.scalar_init)
                {
                    state_ptr(iv, c_phi) += scalar_ptr(i, j, k, 0);
                    state_ptr(iv, c_Pi) += scalar_ptr(i, j, k, 1);
                    state_ptr(iv, c_chi) += scalar_ptr(i, j, k, 2);
                    state_ptr(iv, c_K) += scalar_ptr(i, j, k, 3);
                }

                // Set the entire tensor object here
                if(m_params.tensor_init)
                {
                    state_ptr(iv, c_h11) = hij_ptr(i, j, k, lut[0][0]);
                    state_ptr(iv, c_h12) = hij_ptr(i, j, k, lut[0][1]);
                    state_ptr(iv, c_h13) = hij_ptr(i, j, k, lut[0][2]);
                    state_ptr(iv, c_h22) = hij_ptr(i, j, k, lut[1][1]);
                    state_ptr(iv, c_h23) = hij_ptr(i, j, k, lut[1][2]);
                    state_ptr(iv, c_h33) = hij_ptr(i, j, k, lut[2][2]);

                    state_ptr(iv, c_A11) = Aij_ptr(i, j, k, lut[0][0]);
                    state_ptr(iv, c_A12) = Aij_ptr(i, j, k, lut[0][1]);
                    state_ptr(iv, c_A13) = Aij_ptr(i, j, k, lut[0][2]);
                    state_ptr(iv, c_A22) = Aij_ptr(i, j, k, lut[1][1]);
                    state_ptr(iv, c_A23) = Aij_ptr(i, j, k, lut[1][2]);
                    state_ptr(iv, c_A33) = Aij_ptr(i, j, k, lut[2][2]);
                }
            }
        });
    }
}

/****
    Extraction routines
****/

// Calculates and prints the power spectrum
inline void RandomField::print_power_spectrum(cMultiFab &field_array, SmallDataIO &power_spec_file, const int component = 0)
{ 
    // Set up the isotropic k axis bounds
    Real kiso_max = std::sqrt(3.) * N * M_PI / m_params.L;
    Real dkiso = sqrt(3.)*2.*M_PI/m_params.L;
    Real tolerance = 1.e-12;

    // check the stepping along the diagonal is consistent
    if (kiso_max/dkiso - N/2 > tolerance)
    {
        Error("RandomField::print_power_spectrum Isotropic k axis is too large.");
    }
    // check you aren't sampling above the max sampling rate
    else if (m_params.bin_number > kiso_max/dkiso)
    {
        Error("RandomField::print_power_spectrum Bin number is too large.");
    }
    // check your bin number isn't greater than the max resolvable bins
    else if(m_params.bin_number > m_params.N_readin/2)
    {
        Error("RandomField::print_power_spectrum bin number must be less than N/2.");
    }

    // Set up isotropic k axis and PS map
    Real dk_to_bin = (Real)m_params.bin_number/((Real)N/2);
    Real kmag = 0.;
    Vector<Real> kiso(N/2+1, 0.);

    Vector<Real> ps_map(m_params.bin_number+1, 0.);
    Vector<int> kcount(m_params.bin_number+1, 0);

    for (int s=0; s<=N/2; s++) { kiso[s] = s*dkiso; }

    // Loop to bin the power spectrum at each point
    MFIter::allowMultipleMFIters(true); // Needed to pass the map to the ParallelFor loop
    for (MFIter mfi(field_array); mfi.isValid(); ++mfi) 
    {
        // Make a pointer to the mode function at this MF box
        Array4<GpuComplex<Real>> const& field_ptr = field_array.array(mfi);
        const Box& bx = mfi.fabbox();

        amrex::ParallelFor(bx, [=, &ps_map, &kcount] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            // Check to see if you're in a ghost cell
            IntVect iv{i, j, k};
            bool in_ghost_index = is_ghost_index(iv);
            if(!in_ghost_index)
            {
                Real kmag = get_kmag(iv);

                // make sure you're still in the domain
                if(kmag - kiso_max > tolerance) 
                { 
                    Print() << iv << "\n";
                    Print() << kmag << "," << kiso_max << "\n";
                    Error("RandomField::print_power_spectrum Found magnitude larger than (N/2,N/2,N/2)."); 
                }

                // Loop over the isotropic axis
                for (int s=1; s<=N/2; s++) 
                {
                    // If smaller than the smallest bin
                    if(kmag < kiso[0])
                    {
                        Print() << iv << "\n";
                        Error("RandomField::print_power_spectrum kmag below the kiso domain.");
                    }

                    // If you're larger than the largest bin
                    else if(kmag - kiso[N/2] > tolerance)
                    {
                        Print() << iv << "\n";
                        Error("RandomField::print_power_spectrum kmag above the kiso domain.");
                    }

                    // If you're somewhere in the middle
                    else if (kmag < kiso[s] && kmag >= kiso[(s-1)]) 
                    {
                        Real power = (std::pow(field_ptr(i, j, k, component).real(), 2.0) 
                                    + std::pow(field_ptr(i, j, k, component).imag(), 2.0));

                        if (i != 0 && i != N/2) { power *= 2.; }
                        
                        Gpu::Atomic::Add(&kcount[s-1], 1);
                        if(power > tolerance)
                        {
                            Gpu::Atomic::Add(&ps_map[s-1], power);   
                        }

                        break;
                    }

                    // If you're at the largest bin
                    else if (kmag == kiso[N/2])
                    { 
                        Real power = (std::pow(field_ptr(i, j, k, component).real(), 2.0) 
                                    + std::pow(field_ptr(i, j, k, component).imag(), 2.0));

                        if (i != 0 && i != N/2) { power *= 2.; }
                        
                        Gpu::Atomic::Add(&kcount[N/2], 1);
                        if(power > tolerance)
                        {
                            Gpu::Atomic::Add(&ps_map[N/2], power);
                        }

                        break;
                    }

                    // If you've reached the largest bin but not been captured
                    else if(s > N/2)
                    { 
                        Print() << iv << "\n";
                        Print() << kmag << "\n";
                        Print() << kiso[s] << "," << kiso[s-1] << "\n";
                        Error("RandomField::print_power_spectrum Part of the spectrum isn't captured.");
                    }

                    // If you haven't found the right bin yet
                    else { continue; }
                }
            }
        });
    }

    ParallelAllReduce::Sum(kcount.data(), static_cast<int>(kcount.size()), ParallelContext::CommunicatorSub());
    ParallelAllReduce::Sum(ps_map.data(), static_cast<int>(ps_map.size()), ParallelContext::CommunicatorSub());

    // Print the power spectrum to a new file in data/
#pragma omp single
    for(int s=0; s<=N/2; s++)
    {
        power_spec_file.write_data_line({(kiso[s]+kiso[s+1])/2., (Real)ps_map[s]/kcount[s]});
    }
}

// Finds statistical moment x of given MultiFab
inline Real RandomField::find_field_moment_x(MultiFab &field, const Vector<Real> mean, 
                                             const int moment, const int component)
{
    Real sum = 0.;
    const Real vol = std::pow(N, 3.);

    // Loop over the box to calculate moment x
    for (MFIter mfi(field); mfi.isValid(); ++mfi) 
    {
        Array4<Real> const& field_ptr = field.array(mfi);
        const Box& bx = mfi.fabbox();

        amrex::ParallelFor(bx, [=, &sum] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            sum += std::pow(field_ptr(i, j, k, component) - mean[component], moment);
        });
    }

    ParallelAllReduce::Sum(sum, ParallelContext::CommunicatorSub());
    //if (moment == 3) { Print() << "Components of skewness: ";
    //                   Print() << sum << ", " << sum/vol << "\n"; }

    // Normalise and return moment x
    if (sum == 0) { return 0; }
    else if(moment == 2) { return sqrt(sum/vol); }
    else { return sum/vol; }
}

// Calculates and prints requested moments (any between 1 and 4)
inline Vector<Real> RandomField::print_moment(MultiFab &field, const Vector<std::string> names,  
                                             const Vector<int> &moment_orders, SmallDataIO &file, 
                                             const int is_first_step)
{
    // Trap instance where the user requests too large a moment
    for(const auto moment : moment_orders)
    {
        if(moment > 4) 
        { 
            Error("RandomField::print_moment Chosen moment order has not been implemented");
        }
    }

    // Allocate arrays to store each moment
    const int nc = field.nComp();
    const Real vol = std::pow(N, 3.);
    Vector<Real> means(nc, 0.);
    Vector<Real> stdev(nc, 0.);
    Vector<Real> skew(nc, 0.);
    Vector<Real> kurt(nc, 0.);

    // Find iterators, which determine which moments are requested and their ordering
    Vector<int>::const_iterator start = moment_orders.begin();
    Vector<int>::const_iterator mean_itr = std::find(moment_orders.begin(), moment_orders.end(), 1);
    Vector<int>::const_iterator stdev_itr = std::find(moment_orders.begin(), moment_orders.end(), 2);
    Vector<int>::const_iterator skew_itr = std::find(moment_orders.begin(), moment_orders.end(), 3);
    Vector<int>::const_iterator kurt_itr = std::find(moment_orders.begin(), moment_orders.end(), 4);

    // Allocate vectors to store header line and data lines
    Vector<Real> data_to_print(nc * moment_orders.size(), 0.);
    Vector<std::string> headers(nc * moment_orders.size(), "");

    for (int comp = 0; comp < nc; comp++)
    {
        means[comp] = field.sum(comp)/vol;
        if(mean_itr != moment_orders.end())
        {
            assign_statistics_data(headers, names[comp]+"-mean", data_to_print, means, comp, nc,
                                    mean_itr, start, is_first_step);
        }

        if(moment_orders.back() != 1)
        {
            stdev[comp] = find_field_moment_x(field, means, 2, comp);
            if(stdev_itr != moment_orders.end())
            {
                assign_statistics_data(headers, names[comp]+"-stdev", data_to_print, stdev, comp, nc, 
                                        stdev_itr, start, is_first_step);
            }

            if(moment_orders.back() != 2)
            {
                skew[comp] = find_field_moment_x(field, means, 3, comp);
                skew[comp] /= std::pow(stdev[comp], 3.);

                if(skew_itr != moment_orders.end())
                {
                    assign_statistics_data(headers, names[comp]+"-skew", data_to_print, skew, comp, nc,
                                            skew_itr, start, is_first_step);
                }

                if(moment_orders.back() != 3)
                {
                    kurt[comp] = find_field_moment_x(field, means, 4, comp);
                    kurt[comp] /= std::pow(stdev[comp], 4.);

                    assign_statistics_data(headers, names[comp]+"-kurt", data_to_print, kurt, comp, nc,
                                            kurt_itr, start, is_first_step);
                }
            }
        }
    }

    // Write header and data lines
#pragma omp single
    if(is_first_step) 
    { 
        file.write_header_line(headers); 
    }

#pragma omp single
    file.write_time_data_line(data_to_print);
    
    return stdev;
}

inline void RandomField::derive(const MultiFab &state, MultiFab &out, int dcomp)
{
    BL_PROFILE("RandomField::derive");

    // Extract MultiFab ingredients from state
    BoxArray sba = state.boxArray();
    DistributionMapping sdm = state.DistributionMap();
    MultiFab gij_x(sba, sdm, 6, 0);

    // 0: scalar field
    // 1: conformal factor
    MultiFab scalars_x(sba, sdm, 2, 0);

    // Copy the spatial metric from the state
    Copy(gij_x, state, c_h11, lut[0][0], 1, 0);
    Copy(gij_x, state, c_h12, lut[0][1], 1, 0);
    Copy(gij_x, state, c_h13, lut[0][2], 1, 0);
    Copy(gij_x, state, c_h22, lut[1][1], 1, 0);
    Copy(gij_x, state, c_h23, lut[1][2], 1, 0);
    Copy(gij_x, state, c_h33, lut[2][2], 1, 0);

    int m_c_phi = 0;
    int m_c_chi = 1;
    Copy(scalars_x, state, c_phi, m_c_phi, 1, 0);
    Copy(scalars_x, state, c_chi, m_c_chi, 1, 0);

    // Find background quantities needed to extract \cal R
    const int vol = std::pow(m_params.N_readin, 3);
    const Real K_bar = state.sum(c_K)/vol;
    const Real alpha_bar = state.sum(c_lapse)/vol;
    const Real Pi_bar = state.sum(c_Pi)/vol;
    const Real phi_bar = state.sum(c_phi)/vol;
    const Real chi_bar = state.sum(c_chi)/vol;

    // Remove background from scalar field
    scalars_x.plus(-phi_bar, m_c_phi, 1);
    scalars_x.plus(-chi_bar, m_c_chi, 1);
    scalars_x.mult(1./norm);

    // Undo the normalisation and BSSN-CPT conversion
    for (int l=0; l<3; l++) { gij_x.plus(-1., lut[l][l], 1); }
    gij_x.mult(1./norm);

    // Set up the problem domain in Fourier space
    // And impose that MPI ranks only slice along the i index (for Nyquist conditions)
    IntVect domain_low(0, 0, 0);
    IntVect k_domain_high(N/2, N-1, N-1);
    Box k_domain(domain_low, k_domain_high);
    Array< bool, AMREX_SPACEDIM > const &slicing{true, false, false};
    BoxArray kba = decompose(k_domain, ParallelContext::NProcsAll(), slicing);
    DistributionMapping kdm(kba);

    // Set up the arrays to store the Fourier data sets
    cMultiFab hs_k(kba, kdm, 2, 0);
    cMultiFab gij_k(kba, kdm, 6, 0);
    cMultiFab scalars_k(kba, kdm, 2, 0);
    cMultiFab R_k(kba, kdm, 1, 0);

    hs_k.setVal(0.0);
    gij_k.setVal(0.0);
    scalars_k.setVal(0.0);
    R_k.setVal(0.0);

    // Set up the FFT
    IntVect x_domain_high(N-1, N-1, N-1);
    Box x_domain(domain_low, x_domain_high);
    FFT::R2C<Real> tensor_fft(x_domain, FFT::Info().setBatchSize(gij_k.nComp()));
    FFT::R2C<Real> scalar_fft(x_domain, FFT::Info().setBatchSize(scalars_k.nComp()));

    // Perform the fft
    tensor_fft.forward(gij_x, gij_k);
    scalar_fft.forward(scalars_x, scalars_k);

    // Normalise the fft (fftw style)
    for(int comp = 0; comp < 6; comp++) { gij_k.mult(std::pow(N, -3./2.), comp, 1); }
    for(int comp = 0; comp < 2; comp++) { scalars_k.mult(std::pow(N, -3./2.), comp, 1); }

    // Loop to extract the Fourier-space mode functions
    for (MFIter mfi(gij_k); mfi.isValid(); ++mfi) 
    {
        const Box& bx = mfi.fabbox();
        Array4<GpuComplex<Real>> const& hs_ptr = hs_k.array(mfi);
        Array4<GpuComplex<Real>> const& hij_ptr = gij_k.array(mfi);
        Array4<GpuComplex<Real>> const& scalars_ptr = scalars_k.array(mfi);
        Array4<GpuComplex<Real>> const& R_k_ptr = R_k.array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            IntVect iv{i, j, k};

            if (iv != IntVect{0, 0, 0})
            {
                Vector<Real> mhat = calculate_basis_vector(iv, 0);
                Vector<Real> nhat = calculate_basis_vector(iv, 1);
                Tensor<2, Real> eplus, ecross;

                // Find basis tensors and do the Fourier trick
                for (int l=0; l<3; l++) for (int p=0; p<3; p++)
                {
                    eplus[l][p] = mhat[l]*mhat[p] - nhat[l]*nhat[p];
                    ecross[l][p] = mhat[l]*nhat[p] + nhat[l]*mhat[p];

                    hs_ptr(i, j, k, 0) += (hij_ptr(i, j, k, lut[l][p]) * eplus[l][p])/2.;
                    hs_ptr(i, j, k, 1) += (hij_ptr(i, j, k, lut[l][p]) * ecross[l][p])/2.;
                }

                if (m_params.alpha != 0) { Test_polarisation_tensor_orthonorm(iv, eplus, ecross); }

                // Calculate the TT and scalar-(vector) components of the 
                // metric, by reconstructing hij and subtracting it from \tilde{gamma}_ij
                Tensor<2, GpuComplex<Real>> hij, hSV;
                GpuComplex<Real> hij_tr = 0.;
                GpuComplex<Real> hSV_tr = 0.;
                for (int l=0; l<3; l++) for (int p=0; p<3; p++)
                {
                    hij[l][p] = hs_ptr(i, j, k, 0) * eplus[l][p] + hs_ptr(i, j, k, 1) * ecross[l][p];
                    hSV[l][p] = hij_ptr(i, j, k, lut[l][p]) - hij[l][p];
                }

                // Extract R according to the scheme detailed in 
                // Appendix B (Eq. B1) of arxiv:2502.06783, using hSV as the 
                // spatial metric instead of \tilde{gamma}_ij
                if(m_params.scalar_init)
                {
                    // Find the unitful k vector
                    Vector<Real> iv_k(iv.begin(), iv.end());
                    for(auto& k_comp : iv_k) { k_comp *= 2. * M_PI / m_params.L; }
                    Real kmag = get_kmag(iv);
                    GpuComplex<Real> Phi = 0;

                    // Set the zero mode
                    if(kmag == 0)
                    {
                        R_k_ptr(i, j, k) = GpuComplex<Real>{0., 0.};
                    }

                    else
                    {
                        // converstion from chi and gamma_ij -> Phi
                        for(int l=0; l<3; l++) for(int p=0; p<3; p++)
                        {
                            Phi += (iv_k[l] * iv_k[p] * hSV[l][p])/std::pow(kmag, 2.);
                        }
                        Phi *= 1./4.;
                        Phi += 0.5 * (scalars_ptr(i, j, k, m_c_chi));

                        // Combine the above to find R(k)
                        R_k_ptr(i, j, k) = Phi - K_bar * scalars_ptr(i, j, k, m_c_phi) / alpha_bar / Pi_bar;

                        // Print() << Phi << "\n";
                        // Print() << K_bar << "\n";
                        // Print() << alpha_bar << "\n";
                        // Print() << Pi_bar << "\n";
                        // Print() << scalars_ptr(i, j, k, m_c_phi) << "\n";
                        // Print() << R_k_ptr(i, j, k) << "\n";
                        // Error();
                    }
                }
            }
        });
    }

    apply_nyquist_conditions(hs_k);
    apply_nyquist_conditions(R_k);

    // Make a multifab to store config space mode functions
    // Need to use out to make these ingredients??
    BoxArray xba = out.boxArray();//(x_domain); //
    DistributionMapping xdm = out.DistributionMap();//(xba); //
    MultiFab hs_x(xba, xdm, 2, 0);
    MultiFab R_x(xba, xdm, 1, 0);
    hs_x.setVal(0.0);
    R_x.setVal(0.0);

    // Fourier transform
    FFT::R2C<Real> mode_function_fft(x_domain, FFT::Info().setBatchSize(hs_k.nComp()));
    FFT::R2C<Real> R_fft(x_domain, FFT::Info().setBatchSize(R_x.nComp()));
    mode_function_fft.backward(hs_k, hs_x);
    R_fft.backward(R_k, R_x);

    Test_Parsevals_thm(hs_x, hs_k);
    Test_Parsevals_thm(R_x, R_k);

    // Apply physical normalisation
    hs_x.mult(norm * std::pow(N, -3./2.));
    R_x.mult(norm * std::pow(N, -3./2.));

    for (MFIter mfi(hs_x); mfi.isValid(); ++mfi) 
    {
        Array4<Real> const& out_ptr = out.array(mfi);
        Array4<Real> const& hx_ptr = hs_x.array(mfi);
        Array4<Real> const& Rx_ptr = R_x.array(mfi);

        const Box& bx = mfi.fabbox();
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            const IntVect iv{i, j, k};
            bool in_ghost_index = is_ghost_index(iv);
            if(!in_ghost_index)
            {
                out_ptr(iv, dcomp) = hx_ptr(i, j, k, 0);
                out_ptr(iv, dcomp + 1) = hx_ptr(i, j, k, 1);
                out_ptr(iv, dcomp + 2) = Rx_ptr(i, j, k);
            }
        });
    }
}

// Main extraction routine
inline void RandomField::extract(const MultiFab &state, const std::string data_path, const Real dt,  
                                 const Real cur_time, const int restart_time, const int first_step)
{
    BL_PROFILE("RandomField::extract");

    // Extract MultiFab ingredients from state
    BoxArray sba = state.boxArray();
    DistributionMapping sdm = state.DistributionMap();
    MultiFab gij_x(sba, sdm, 6, 0);

    // 0: scalar field
    // 1: conformal factor
    MultiFab scalars_x(sba, sdm, 2, 0);

    // Copy the spatial metric from the state
    Copy(gij_x, state, c_h11, lut[0][0], 1, 0);
    Copy(gij_x, state, c_h12, lut[0][1], 1, 0);
    Copy(gij_x, state, c_h13, lut[0][2], 1, 0);
    Copy(gij_x, state, c_h22, lut[1][1], 1, 0);
    Copy(gij_x, state, c_h23, lut[1][2], 1, 0);
    Copy(gij_x, state, c_h33, lut[2][2], 1, 0);

    int m_c_phi = 0;
    int m_c_chi = 1;
    Copy(scalars_x, state, c_phi, m_c_phi, 1, 0);
    Copy(scalars_x, state, c_chi, m_c_chi, 1, 0);

    // Find background quantities needed to extract \cal R
    const int vol = std::pow(m_params.N_readin, 3);
    const Real K_bar = state.sum(c_K)/vol;
    const Real alpha_bar = state.sum(c_lapse)/vol;
    const Real Pi_bar = state.sum(c_Pi)/vol;
    const Real phi_bar = state.sum(c_phi)/vol;
    const Real chi_bar = state.sum(c_chi)/vol;

    // Remove background from scalar field
    scalars_x.plus(-phi_bar, m_c_phi, 1);
    scalars_x.plus(-chi_bar, m_c_chi, 1);
    scalars_x.mult(1./norm);

    // Undo the normalisation and BSSN-CPT conversion
    for (int l=0; l<3; l++) { gij_x.plus(-1., lut[l][l], 1); }
    gij_x.mult(1./norm);

    // Set up the problem domain in Fourier space
    // And impose that MPI ranks only slice along the i index (for Nyquist conditions)
    IntVect domain_low(0, 0, 0);
    IntVect k_domain_high(N/2, N-1, N-1);
    Box k_domain(domain_low, k_domain_high);
    Array< bool, AMREX_SPACEDIM > const &slicing{true, false, false};
    BoxArray kba = decompose(k_domain, ParallelContext::NProcsAll(), slicing);
    DistributionMapping kdm(kba);

    // Set up the arrays to store the Fourier data sets
    cMultiFab hs_k(kba, kdm, 2, 0);
    cMultiFab gij_k(kba, kdm, 6, 0);
    cMultiFab scalars_k(kba, kdm, 2, 0);
    cMultiFab R_k(kba, kdm, 1, 0);

    hs_k.setVal(0.0);
    gij_k.setVal(0.0);
    scalars_k.setVal(0.0);
    R_k.setVal(0.0);

    // Set up the FFT
    IntVect x_domain_high(N-1, N-1, N-1);
    Box x_domain(domain_low, x_domain_high);
    FFT::R2C<Real> tensor_fft(x_domain, FFT::Info().setBatchSize(gij_k.nComp()));
    FFT::R2C<Real> scalar_fft(x_domain, FFT::Info().setBatchSize(scalars_k.nComp()));

    // Perform the fft
    tensor_fft.forward(gij_x, gij_k);
    scalar_fft.forward(scalars_x, scalars_k);

    // Normalise the fft (fftw style)
    for(int comp = 0; comp < 6; comp++) { gij_k.mult(std::pow(N, -3./2.), comp, 1); }
    for(int comp = 0; comp < 2; comp++) { scalars_k.mult(std::pow(N, -3./2.), comp, 1); }

    // Set variables to store the maximum trace 
    // of the scalar and tensor components
    Real hij_tr_max = 0.;
    Real hSV_tr_max = 0.;

    // Loop to extract the Fourier-space mode functions
    for (MFIter mfi(gij_k); mfi.isValid(); ++mfi) 
    {
        const Box& bx = mfi.fabbox();
        Array4<GpuComplex<Real>> const& hs_ptr = hs_k.array(mfi);
        Array4<GpuComplex<Real>> const& hij_ptr = gij_k.array(mfi);
        Array4<GpuComplex<Real>> const& scalars_ptr = scalars_k.array(mfi);
        Array4<GpuComplex<Real>> const& R_k_ptr = R_k.array(mfi);

        amrex::ParallelFor(bx, [=, &hij_tr_max, &hSV_tr_max] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            IntVect iv{i, j, k};
            
            Vector<Real> mhat = calculate_basis_vector(iv, 0);
            Vector<Real> nhat = calculate_basis_vector(iv, 1);
            Tensor<2, Real> eplus, ecross;

            // Find basis tensors and do the Fourier trick
            for (int l=0; l<3; l++) for (int p=0; p<3; p++)
            {
                eplus[l][p] = mhat[l]*mhat[p] - nhat[l]*nhat[p];
                ecross[l][p] = mhat[l]*nhat[p] + nhat[l]*mhat[p];

                hs_ptr(i, j, k, 0) += (hij_ptr(i, j, k, lut[l][p]) * eplus[l][p])/2.;
                hs_ptr(i, j, k, 1) += (hij_ptr(i, j, k, lut[l][p]) * ecross[l][p])/2.;
            }

            if (m_params.alpha != 0) { Test_polarisation_tensor_orthonorm(iv, eplus, ecross); }

            // Calculate the TT and scalar-(vector) components of the 
            // metric, by reconstructing hij and subtracting it from \tilde{gamma}_ij
            Tensor<2, GpuComplex<Real>> hij, hSV;
            GpuComplex<Real> hij_tr = 0.;
            GpuComplex<Real> hSV_tr = 0.;
            for (int l=0; l<3; l++) for (int p=0; p<3; p++)
            {
                hij[l][p] = hs_ptr(i, j, k, 0) * eplus[l][p] + hs_ptr(i, j, k, 1) * ecross[l][p];
                hSV[l][p] = hij_ptr(i, j, k, lut[l][p]) - hij[l][p];

                // Find the traces of these components, as a diagnostic
                if(l==p) { hij_tr += hij[l][p]; hSV_tr += hSV[l][p]; }
            }
            
            Real hij_tr_mag = sqrt(pow(hij_tr.real(), 2.) + pow(hij_tr.imag(), 2.));
            Real hSV_tr_mag = sqrt(pow(hSV_tr.real(), 2.) + pow(hSV_tr.imag(), 2.));

            if (hij_tr_mag > hij_tr_max) { hij_tr_max = hij_tr_mag; }
            if (hSV_tr_mag > hSV_tr_max) { hSV_tr_max = hSV_tr_mag; }

            // Confirm hij is trace-free in Fourier space
            if (std::abs(hij_tr_mag) > tolerance)
            {
                Print() << iv << ": " << hij_tr_mag << "\n";
                Error("RandomField::extract, "
                      "hij trace magnitude too large in extraction");
            }

            // Extract R according to the scheme detailed in 
            // Appendix B (Eq. B1) of arxiv:2502.06783, using hSV as the 
            // spatial metric instead of \tilde{gamma}_ij
            if(m_params.scalar_init)
            {
                // Find the unitful k vector
                Vector<Real> iv_k(iv.begin(), iv.end());
                for(auto& k_comp : iv_k) { k_comp *= 2. * M_PI / m_params.L; }
                Real kmag = get_kmag(iv);
                GpuComplex<Real> Phi = 0;

                // Set the zero mode
                if(kmag == 0)
                {
                    R_k_ptr(i, j, k) = GpuComplex<Real>{0., 0.};
                }

                else
                {
                    // converstion from chi and gamma_ij -> Phi
                    for(int l=0; l<3; l++) for(int p=0; p<3; p++)
                    {
                        Phi += (iv_k[l] * iv_k[p] * hSV[l][p])/std::pow(kmag, 2.);
                    }
                    Phi *= 1./4.;
                    Phi += 0.5 * (scalars_ptr(i, j, k, m_c_chi));

                    // Combine the above to find R(k)
                    R_k_ptr(i, j, k) = Phi - K_bar * scalars_ptr(i, j, k, m_c_phi) / alpha_bar / Pi_bar;
                }
            }
        });
    }

    // Output the max traces of the tensor components as a diagnostic
    SmallDataIO trace_file(data_path+"tensor-traces", dt, cur_time, 
                                restart_time, SmallDataIO::APPEND, first_step, ".dat");
    if(first_step) 
    { 
        trace_file.write_header_line({"hij trace max", "hSV trace max"}); 
    }
    trace_file.write_time_data_line({hij_tr_max, hSV_tr_max}); 

    apply_nyquist_conditions(hs_k);
    apply_nyquist_conditions(R_k);

    // Find the binned PS for each mode function and print to data/
    if((m_params.calc_binned_power_spectrum) 
	    && (static_cast<int>(cur_time/dt) % m_params.plot_int == 0))
    {
	    Print() << "RandomField::extract, Time step at print: ";
        Print() << static_cast<int>(std::round(cur_time/dt)) << "\n";

        std::string spec_path = make_subdirectory(data_path, "spectra", first_step);
        Vector<std::string> filenames(2, "");

        for(int comp = 0; comp < hs_k.nComp(); comp++)
        {
            filenames[comp] = spec_path+"spectrum-comp-"+std::to_string(comp)+"-time-";
            SmallDataIO spectrum_file(filenames[comp], dt, cur_time, restart_time, SmallDataIO::NEW, first_step, ".dat");
            print_power_spectrum(hs_k, spectrum_file, comp);
        }

        std::string filename = spec_path+"spectrum-Rk-time-";
        SmallDataIO spectrum_file(filename, dt, cur_time, restart_time, SmallDataIO::NEW, first_step, ".dat");
        print_power_spectrum(R_k, spectrum_file, 0);
    }

    // Find mode functions in configuration space if requested,
    // and find the statistics (orders 1-4) of the polarisation fields 
    // and the R field. 
    // And print the tensor-to-scalar ratio if requested.
    if(m_params.calc_higher_order_statistics)
    {
        // Make a multifab to store config space mode functions
        BoxArray xba(x_domain);
        DistributionMapping xdm(xba);
        MultiFab hs_x(xba, xdm, 2, 0);
        MultiFab R_x(xba, xdm, 1, 0);
        hs_x.setVal(0.0);
        R_x.setVal(0.0);

        // Fourier transform
        FFT::R2C<Real> tensor_mode_function_fft(x_domain, FFT::Info().setBatchSize(hs_k.nComp()));
        FFT::R2C<Real> scalar_mode_function_fft(x_domain, FFT::Info().setBatchSize(R_k.nComp()));
        tensor_mode_function_fft.backward(hs_k, hs_x);
        scalar_mode_function_fft.backward(R_k, R_x);

        if (m_params.tensor_init) 
        { 
            for (int c=0; c<hs_x.nComp(); c++)
            {
                MultiFab hx_tmp(hs_x, amrex::make_alias, c, 1);
                cMultiFab hk_tmp(hs_k, amrex::make_alias, c, 1);
                Test_Parsevals_thm(hx_tmp, hk_tmp);
            } 
        }
        if (m_params.scalar_init) { Test_Parsevals_thm(R_x, R_k); }

        // Apply physical normalisation
        hs_x.mult(norm * std::pow(N, -3./2.));
        R_x.mult(norm * std::pow(N, -3./2.));

        // Print() << "Max tensor polarisations: " << hs_x.max(0) << ", " << hs_x.max(1) << "\n";
        // Print() << "R max and min bounds: " << R_x.max(0) << ", " << R_x.min(0) << "\n";

        int output_comps = hs_x.nComp() + R_x.nComp();
        MultiFab out_MF(hs_x.boxArray(), hs_x.DistributionMap(), output_comps, 0);
        Copy(out_MF, R_x, 0, 0, R_k.nComp(), 0);
        Copy(out_MF, hs_x, 0, R_k.nComp(), hs_x.nComp(), 0);

        // Print mode functions if requested
        if(m_params.print_mode_functions == 1)
        {
            std::string mf_path = make_subdirectory(data_path, "mode-functions", first_step);
            std::string filename = mf_path+"mode-function-"+std::to_string(cur_time/dt);

            for (MFIter mfi(hs_x); mfi.isValid(); ++mfi) 
            {
                Array4<Real> const& hx_ptr = hs_x.array(mfi);
                const Box& bx = mfi.fabbox();

                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    AllPrintToFile(filename) << i + N*(j + N*k) << ",";
                    for(int c=0; c<2; c++)
                    {
                        AllPrintToFile(filename).SetPrecision(14) << hx_ptr(i, j, k, c) << ",";
                    }
                    AllPrintToFile(filename) << "\n";
                });
            }
        }

        // Calculate and print field moments if requested
        Vector<Real> stdevs;
        if (m_params.calc_higher_order_statistics)
        {
            SmallDataIO stats_file(data_path+"field-statistics", dt, cur_time, 
                                    restart_time, SmallDataIO::APPEND, first_step, ".dat");

            if (!m_params.orders.empty())
            {
                Vector<std::string> names{"R","hplus","hcross"};
                stdevs = print_moment(out_MF, names, m_params.orders, stats_file, first_step);
                // Print() << stdevs[0] << "\n";
                // Error();
            }
        }
        
        SmallDataIO ts_file(data_path+"tensor-scalar-ratio", dt, cur_time, 
                                    restart_time, SmallDataIO::APPEND, first_step, ".dat");
#pragma omp single
        if(first_step) 
    	{ 
            ts_file.write_header_line({"T/S ratio (plus)", "T/S ratio (cross)"}); 
        }
        
#pragma omp single
    	ts_file.write_time_data_line({stdevs[1] / stdevs[0], stdevs[2] / stdevs[0]}); 
    }
}

#endif /* RANDOMFIELD_IMPL_HPP_*/

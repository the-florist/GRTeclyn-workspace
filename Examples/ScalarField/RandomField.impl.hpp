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
            Print() << "Directory creation failed for " << new_path << "\n";
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

/****
    Tests
****/

inline void RandomField::Test_is_trace_free(MultiFab &field)
{
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
                Print() << iv << ": " << sum;
                Error("RandomField::Test_is_trace_free Trace-free test failed here.");
            }
        });
    }
    
}

/****
    Initialisation routines
****/

// Generate unique random draws for each MFI box.
inline void RandomField::make_random_draws(auto &rand_fab, Box &domain)
{
    BoxArray ba = rand_fab.boxArray();
    DistributionMapping dm = rand_fab.DistributionMap();
    FabArray<BaseFab<GpuArray<Real, 4>>> tmp(ba, dm, 1, 0, MFInfo{}.SetArena(The_Cpu_Arena()));

    for(MFIter mfi(tmp); mfi.isValid(); ++mfi)
    {
        Box const& bx = mfi.validbox();
        auto const& tmp_ptr = tmp.array(mfi);

        std::mt19937 generator;
        std::uniform_real_distribution<Real> distribution(Real(0), Real(1));

        auto offset = domain.index(bx.smallEnd()) * 4;
        for(int ofs = 0; ofs < offset; ofs++)
        {
            distribution(generator);
        }
        amrex::LoopOnCpu(bx, [&] (int i, int j, int k)
        {
            auto &field_point = tmp_ptr(i, j, k);
            for(int l=0; l<4; l++)
            {
                field_point[l] = distribution(generator);
            }
        });
    }

    rand_fab.ParallelCopy(tmp);
}

// Returns analytic power spectrum in modulus/argument form
inline GpuComplex<Real> RandomField::calculate_mode_function(const double km, const std::string spec_type)
{
    // Deals with k=0 case, which is undefined if m=0
    if(km < 1.e-23) { return 0.; }
    
    // Stores modulus and argument 
    Real ms_mag = 0.;
    Real ms_arg = 0.;

    // Hubble at t=0, needed for tensor solution
    Real H0 = sqrt((4.0 * M_PI/3.0/pow(m_params.Mp, 2.))
                * (pow(m_background_params.m * m_background_params.phi0, 2.0) 
                    + pow(m_background_params.Pi0, 2.)));

    double kpr = km/H0;
    if (spec_type == "position") // Position mode funcion
    {
        ms_mag = sqrt((1.0/km + H0*H0/pow(km, 3.))/2.);
        ms_arg = atan2((cos(kpr) + kpr*sin(kpr)), (kpr*cos(kpr) - sin(kpr)));
    }
    else if (spec_type == "velocity") // Velocity mode funcion
    {
        ms_mag = sqrt(km/2.);
        ms_arg = -atan2(cos(kpr), sin(kpr));
    }
    else { Error("RandomField::calculate_mode_function Value of spec_type not allowed."); }

    // Construct the mode function and return it
    GpuComplex<Real> ps(ms_mag * cos(ms_arg), ms_mag * sin(ms_arg));
    return ps;
}

// Turns analytic PS into GRF and applies window function if requested
inline GpuComplex<Real> RandomField::calculate_random_field(const IntVect iv, const std::string spectrum_type, 
                                                            const Real rand_amp, const Real rand_phase)
{
    GpuComplex<Real> value(0., 0.);

    // Find kmag with FFTW-style inversion on the last two indices
    int i = iv[0];
    int j = invert_index(iv[1]);
    int k = invert_index(iv[2]);

    double kmag = std::sqrt(i*i + j*j + k*k) * 2 * M_PI / L;

    // Find the analytic power spectrum
    value = calculate_mode_function(kmag, spectrum_type);

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
        double ks = m_params.kstar * 2. * M_PI/L;
        double Dt = L/m_params.Delta;
        value *= 0.5 * (1.0 - tanh(Dt * (kmag - ks))); 
    }

    return value;
}

// Calculates basis vectors required for polarisation tensors
inline Vector<Real> RandomField::calculate_basis_vector(const IntVect iv, const int which_vector)
{
    // FFTW-style inversion with sign on the last two indices
    int i = iv[0];
    int j = invert_index_with_sign(iv[1]);
    int k = invert_index_with_sign(iv[2]);

    Vector<Real> mhat(3, 0.);
    Vector<Real> nhat(3, 0.);

    if (i > 0.) 
    {
        if (k == 0. && j == 0.) { mhat[0] = 1.; mhat[1] = 0.; mhat[2] = 0.; 
                                  nhat[0] = 0.; nhat[1] = 1.; nhat[2] = 0.; 
                                }

        else { mhat[0] = j/sqrt(k*k+j*j); mhat[1] = -k/sqrt(k*k+j*j); mhat[2] = 0.L;
               nhat[0] = k*i/sqrt(i*i*(k*k + j*j) + pow(k*k + j*j, 2.));
               nhat[1] = i*j/sqrt(i*i*(k*k + j*j) + pow(k*k + j*j, 2.));
               nhat[2] = -(k*k + j*j)/sqrt(i*i*(k*k + j*j) + pow(k*k + j*j, 2.)); 
             }
    }

    else if (std::abs(j) > 0) { mhat[0] = 0.; mhat[1] = 0.; mhat[2] = -1.;
                      nhat[0] = -j/sqrt(j*j + k*k);
                      nhat[1] = k/sqrt(j*j + k*k);
                      nhat[2] = 0.; 
                    }

    else if (std::abs(k) > 0) { mhat[0] = 0.; mhat[1] = 1.; mhat[2] = 0.;
                      nhat[0] = 0.; nhat[1] = 0.; nhat[2] = 1.;
                    }

    else if (i==0 && j==0 && k==0) { ; }

    else 
    {
        Error("RandomField::calculate_polarisation_tensors Part of Fourier grid not covered.");
    }

    if(which_vector == 0) { return mhat; }
    else if(which_vector == 1) { return nhat; }
    else { Error("RandomField::calculate_basis_vector Incompatable vector type."); return mhat; }
}

// Assembles full tensor initial conditions given two mode functions
inline GpuComplex<Real> RandomField::calculate_tensor_initial_conditions(const IntVect iv, const int l, const int p, 
                                                                         const GpuComplex<Real> plus_field, 
                                                                         const GpuComplex<Real> cross_field)
{
    Vector<Real> mhat(3, 0.);
    Vector<Real> nhat(3, 0.);

    mhat = calculate_basis_vector(iv, 0);
    nhat = calculate_basis_vector(iv, 1);

    // Assemble the polarisation tensors
    Real eplus = mhat[l]*mhat[p] - nhat[l]*nhat[p];
    Real ecross = mhat[l]*nhat[p] + nhat[l]*mhat[p];

    return (eplus * plus_field + ecross * cross_field)/std::sqrt(2.);
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
    InitRandom(m_params.random_seed);

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

    // Make the Fourier transform
    IntVect x_domain_high(N-1, N-1, N-1);
    Box x_domain(domain_low, x_domain_high);
    FFT::R2C<Real> random_field_fft(x_domain, FFT::Info().setBatchSize(hij_k.nComp()));

    FabArray<BaseFab<GpuArray<Real, 4>>> random_draws(kba, kdm, 1, 0);
    make_random_draws(random_draws, k_domain);

    std::string Filename = "/nfs/st01/hpc-gr-epss/eaf49/GRTeclyn-dump/hs-k-init";
    for (MFIter mfi(hs_k); mfi.isValid(); ++mfi) 
    {
        // Define the domain on this MPI rank
        const Box& bx = mfi.fabbox();
        auto const& random_box_ptr = random_draws.const_array(mfi);
        //int count = 0;

        // Make a pointer to the mode functions at this MF box
        Array4<GpuComplex<Real>> const& hs_ptr = hs_k.array(mfi);
        Array4<GpuComplex<Real>> const& hij_ptr = hij_k.array(mfi);

        Array4<GpuComplex<Real>> const& As_ptr = As_k.array(mfi);
        Array4<GpuComplex<Real>> const& Aij_ptr = Aij_k.array(mfi);

        // Loop to create mode functions, then hij(k) and Aij(k)
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            IntVect iv = {i, j, k};
            auto const& random_field_ptr = random_box_ptr(i, j, k);

            // Find the mode function realisation
            for(int p=0; p<2; p++)
            {
                Real draw1 = amrex::Random();//random_field_ptr[2*p];
                Real draw2 = amrex::Random();//random_field_ptr[2*p+1];

                /*if(count==0)
                {
                    AllPrint() << ParallelContext::MyProcSub() << "," << draw1 << "\n";
                    count++;
                }*/

                hs_ptr(i, j, k, p) = calculate_random_field(iv, "position", draw1, draw2);
                As_ptr(i, j, k, p) = calculate_random_field(iv, "velocity", draw1, draw2);
            }

            // Find basis tensors and initial tensor realisation
            for (int l=0; l<3; l++) for (int p=0; p<3; p++)
            {
                hij_ptr(i, j, k, lut[l][p]) = calculate_tensor_initial_conditions(iv, l, p, hs_ptr(i, j, k, 0), hs_ptr(i, j, k, 1));
                Aij_ptr(i, j, k, lut[l][p]) = calculate_tensor_initial_conditions(iv, l, p, As_ptr(i, j, k, 0), As_ptr(i, j, k, 1));
            }
        });
    }

    // Apply the Nyquist conditions
    apply_nyquist_conditions(hs_k);
    apply_nyquist_conditions(hij_k);
    apply_nyquist_conditions(Aij_k);

    // Do the Fourier transform
    random_field_fft.backward(hij_k, hij_x);
    random_field_fft.backward(Aij_k, Aij_x);

    // Apply normalisation into physical units
    hij_x.mult(norm);
    Aij_x.mult(norm);

    // Test is trace-free
    Test_is_trace_free(hij_x);
    Test_is_trace_free(Aij_x);

    // Convert to BSSN variables using the BSSN-CPT dictionary
    for (int l=0; l<3; l++) { hij_x.plus(1., lut[l][l], 1); }
    Aij_x.mult(-0.5);

    // Put these initial conditions into state
    for (MFIter mfi(hij_x); mfi.isValid(); ++mfi) 
    {
        Array4<Real> const& state_ptr = state.array(mfi);
        Array4<Real> const& hij_ptr = hij_x.array(mfi);
        Array4<Real> const& Aij_ptr = Aij_x.array(mfi);

        const Box& bx = mfi.fabbox();
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            const IntVect iv{i, j, k};
            bool in_ghost_index = is_ghost_index(iv);
            if(!in_ghost_index)
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
        });
    }
}

#endif /* RANDOMFIELD_IMPL_HPP_*/

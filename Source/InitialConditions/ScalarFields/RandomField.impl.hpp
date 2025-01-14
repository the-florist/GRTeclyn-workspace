/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */


#if !defined(RANDOMFIELD_HPP_)
#error "This file should only be included via RandomField.hpp"
#endif

#ifndef RANDOMFIELD_IMPL_HPP_
#define RANDOMFIELD_IMPL_HPP_

// Used for symmetry condition
inline int RandomField::flip_index(int indx) { return std::abs(N - indx); }

// Used for symmetry condition and calculation of kmag
inline int RandomField::invert_index(int indx) { return (int)(N/2 - std::abs(N/2 - indx)); }

// Used for calculation of polarisation tensors
inline int RandomField::invert_index_with_sign(int indx) 
{ 
    if(indx <= N/2) { return indx; }
    else { return std::abs(N/2 - indx) - N/2; }
}

/****
    Initialisation routines
****/

inline GpuComplex<Real> RandomField::calculate_mode_function(double km, std::string spec_type)
{
    // Deals with k=0 case, undefined if m=0
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
    else { Error("RandomField::calculate_power_spectrum Value of spec_type not allowed."); }

    // Construct the mode function and return it
    GpuComplex<Real> ps(ms_mag * cos(ms_arg), ms_mag * sin(ms_arg));
    return ps;
}

inline GpuComplex<Real> RandomField::calculate_random_field(int i, int J, int K, std::string spectrum_type, 
                                                                Real rand_amp, Real rand_phase)
{
    // Storage for the returned value
    GpuComplex<Real> value(0., 0.);

    // Find kmag with FFTW-style inversion on the first two indices
    int j = invert_index(J);
    int k = invert_index(K);
    double kmag = std::sqrt(i*i + j*j + k*k) * 2 * M_PI / m_params.L;

    // Find the linearised solution
    value = calculate_mode_function(kmag, spectrum_type);

    // Add stochastic perturbations
    if(m_params.use_rand == 1)
    {
        BL_PROFILE("RandomField::calculate_random_field Random initialisation is used");

        // Make one random draw for the amplitude and phase
        Real rand_mod = sqrt(-2. * log(rand_amp)); // Rayleigh distribution about |h|
        Real rand_arg = 2. * M_PI * rand_phase;      // Uniform random phase

        // Multiply amplitude by Rayleigh draw
        value *= rand_mod;

        // Apply the random phase, assuming MS phase is accounted for
        Real new_real = value.real() * cos(rand_arg) - value.imag() * sin(rand_arg);
        Real new_imag = value.real() * sin(rand_arg) + value.imag() * cos(rand_arg);
        GpuComplex<Real> new_value(new_real, new_imag);
	
        value = new_value;
    }

    // Apply a window function
    if(m_params.use_window == 1) 
    { 
        BL_PROFILE("RandomField::calculate_random_field Window function is used")
        double ks = m_params.kstar * 2. * M_PI/m_params.L;
        double Dt = m_params.L/m_params.Delta;
        value *= 0.5 * (1. - tanh(Dt * (kmag - ks))); 
    }

    return value;
}

inline GpuComplex<Real> RandomField::calculate_tensor_initial_conditions(int i, int J, int K, int l, int p, 
                                        GpuComplex<Real> plus_field, GpuComplex<Real> cross_field)
{
    // Find kmag with FFTW-style inversion on the first two indices
    int j = invert_index_with_sign(J);
    int k = invert_index_with_sign(K);

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

    Real eplus = 0.;
    Real ecross = 0.;
    eplus = mhat[l]*mhat[p] - nhat[l]*nhat[p];
    ecross = mhat[l]*nhat[p] + nhat[l]*mhat[p];

    return (eplus * plus_field + ecross * cross_field)/std::sqrt(2.);
}

inline void RandomField::apply_nyquist_conditions(int i, int j, int k, Array4<GpuComplex<Real>> const& field)
{
    // Nyquist node condition
    if ((i==0 || i==N/2) && (j==0 || j==N/2) && (k==0 || k== N/2))
    {
        for(int comp = 0; comp < field.nComp(); comp++)
        {
            GpuComplex<Real> temp(field(i, j, k, comp).real(), 0.);
            field(i, j, k, comp) = temp;
        }
    }

    // Nyquist axis condition
    if (i==0 || i==N/2) 
    {
        if((k>N/2 && j==N/2) || (k==0 && j>N/2) || (k>N/2 && j==0) || (k==N/2 && j>N/2))
        {
            for(int comp = 0; comp < field.nComp(); comp++) 
            {
                GpuComplex<Real> temp(field(i, invert_index(j), invert_index(k), comp).real(), 
                                        -field(i, invert_index(j), invert_index(k), comp).imag());
                field(i, j, k, comp) = temp;
            }
        }
        else if(j > N/2)
        {
            for(int comp = 0; comp < field.nComp(); comp++) 
            {
                GpuComplex<Real> temp(field(i, invert_index(j), flip_index(k), comp).real(), 
                                        -field(i, invert_index(j), flip_index(k), comp).imag());
                field(i, j, k, comp) = temp;
            }
        }
    }
}

inline void RandomField::assign_to_grid(const CellData<Real> &current_cell, const Real tensor_field) const
{
    current_cell[c_h11] = tensor_field;
}

inline void RandomField::init(amrex::MultiFab &state)
{
    BL_PROFILE("RandomField::init_random_field");
    InitRandom(m_params.random_seed);

    // Set up the problem domain and MF ingredients (Real space)
    IntVect domain_low(0, 0, 0);
    IntVect domain_high(N-1, N-1, N-1);
    Box domain(domain_low, domain_high);
    BoxArray xba(domain);
    DistributionMapping xdm(xba);

    // Make the fft and store the problem domain and MF ingredients (Fourier space)
    FFT::R2C<Real> random_field_fft(domain);
    auto const& [kba, kdm] = random_field_fft.getSpectralDataLayout();

    // Set up the arrays to store the in/out data sets
    cMultiFab hs_k(kba, kdm, 2, 0);
    cMultiFab As_k(kba, kdm, 2, 0);
    //MultiFab hs_x(xba, xdm, 2, 0);

    cMultiFab hij_k(kba, kdm, 6, 0);
    MultiFab hij_x(xba, xdm, 6, 0);
    cMultiFab Aij_k(kba, kdm, 6, 0);
    MultiFab Aij_x(xba, xdm, 6, 0);

    std::string Filename = "/nfs/st01/hpc-gr-epss/eaf49/GRTeclyn-dump/GRTeclyn-hij-k";
    
    // Loop to create Fourier-space tensor object
    for (MFIter mfi(hs_k); mfi.isValid(); ++mfi) 
    {
        // Make a pointer to the mode functions at this MF box
        Array4<GpuComplex<Real>> const& hs_ptr = hs_k.array(mfi);
        Array4<GpuComplex<Real>> const& hij_ptr = hij_k.array(mfi);

        Array4<GpuComplex<Real>> const& As_ptr = As_k.array(mfi);
        Array4<GpuComplex<Real>> const& Aij_ptr = Aij_k.array(mfi);

        const Box& bx = mfi.fabbox();

        // Loop to create mode functions then hij(k)
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            // Find the mode function realisation
            for(int p=0; p<2; p++)
            {
                Real draw1 = amrex::Random();
                Real draw2 = amrex::Random();

                hs_ptr(i, j, k, p) = calculate_random_field(i, j, k, "position", draw1, draw2);
                As_ptr(i, j, k, p) = calculate_random_field(i, j, k, "velocity", draw1, draw2);
            }

            // Find basis tensors and initial tensor realisation
            for (int l=0; l<3; l++) for (int p=l; p<3; p++)
            {
                hij_ptr(i, j, k, lut[l][p]) = calculate_tensor_initial_conditions(i, j, k, l, p, 
                                                hs_ptr(i, j, k, 0), hs_ptr(i, j, k, 1));
                Aij_ptr(i, j, k, lut[l][p]) = calculate_tensor_initial_conditions(i, j, k, l, p, 
                                                As_ptr(i, j, k, 0), As_ptr(i, j, k, 1));
            }
        });

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            //apply_nyquist_conditions(i, j, k, hs_ptr);
            apply_nyquist_conditions(i, j, k, hij_ptr);
            apply_nyquist_conditions(i, j, k, Aij_ptr);
        });
    }

    for(int fcomp = 0; fcomp < hij_k.nComp(); fcomp++)
    {
        cMultiFab hij_k_slice(hij_k, make_alias, fcomp, 1);
        MultiFab hij_x_slice(hij_x, make_alias, fcomp, 1);
        random_field_fft.backward(hij_k_slice, hij_x_slice);

        cMultiFab Aij_k_slice(Aij_k, make_alias, fcomp, 1);
        MultiFab Aij_x_slice(Aij_x, make_alias, fcomp, 1);
        random_field_fft.backward(Aij_k_slice, Aij_x_slice);
    }

    hij_x.mult(norm);
    Aij_x.mult(norm);

    for (int l=0; l<3; l++) { hij_x.plus(1., lut[l][l], 1); }
    Aij_x.mult(-0.5);

    auto const &state_array = state.arrays();

    amrex::ParallelFor(
        state, state.nGrowVect(),
        [=] AMREX_GPU_DEVICE(int box_ind, int i, int j, int k) noexcept
        {
            const IntVect iv{i, j, k};
            const Array4<Real> a_array = state_array[box_ind];
            a_array(iv, c_h11) = 0.;
        });

    //Add(state, hij_x, c_h11, lut[0][0], 1, 0);

    //print_tensor_moment(1, hij_x);

    //hx = &hij_x;
    //(*hx).nComp();

    /*for (MFIter mfi(hs_k); mfi.isValid(); ++mfi) 
    {
        Array4<Real> const& hx_ptr = (*hx).array(mfi);
        const Box& bx = mfi.fabbox();
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            std::cout << "Inside init() now...\n";
            std::cout << hx_ptr(i, j, k, 0) << "\n";
            Error();
        });
    }*/
}

/*template <class data_t>
void RandomField::compute(int i, int j, int k,
                           const amrex::Array4<data_t> &state) const
{
    std::cout << "Inside compute now...\n";
    std::cout << *hx(i, j, k, 0) << "\n";
    Error();
}*/

/****
    Extraction routines
****/

inline void RandomField::print_tensor_moment(int moment_order, MultiFab &field)
{
    std::cout << field.nComp() << "\n";
    Error();
}

#endif /* RANDOMFIELD_IMPL_HPP_*/

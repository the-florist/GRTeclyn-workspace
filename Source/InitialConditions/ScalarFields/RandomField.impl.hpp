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

inline GpuComplex<Real> RandomField::calculate_random_field(int I, int J, int k, std::string spectrum_type)
{
    // Storage for the returned value
    GpuComplex<Real> value(0., 0.);

    // Find kmag with FFTW-style inversion on the first two indices
    int i = invert_index(I);
    int j = invert_index(J);
    double kmag = std::sqrt(i*i + j*j + k*k) * 2 * M_PI / m_params.L;

    // Find the linearised solution
    value = calculate_mode_function(kmag, spectrum_type);

    // Add stochastic perturbations
    if(m_params.use_rand == 1)
    {
        BL_PROFILE("RandomField::calculate_random_field Random initialisation is used")
        // Make one random draw for the amplitude and phase
        Real rand_mod = sqrt(-2. * log(amrex::Random())); // Rayleigh distribution about |h|
        Real rand_arg = 2. * M_PI * amrex::Random();      // Uniform random phase

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

inline Real RandomField::basis_vector(int I, int J, int k, int l, int which)
{
    // Find kmag with FFTW-style inversion on the first two indices
    int i = invert_index_with_sign(I);
    int j = invert_index_with_sign(J);

    Vector<Real> mhat(3, 0.);
    Vector<Real> nhat(3, 0.);

    if (k > 0.) 
    {
        if (i == 0. && j == 0.) { mhat[0] = 1.; mhat[1] = 0.; mhat[2] = 0.; 
                                  nhat[0] = 0.; nhat[1] = 1.; nhat[2] = 0.; 
                                }

        else { mhat[0] = j/sqrt(i*i+j*j); mhat[1] = -i/sqrt(i*i+j*j); mhat[2] = 0.L;
               nhat[0] = k*i/sqrt(k*k*(i*i + j*j) + pow(i*i + j*j, 2.));
               nhat[1] = k*j/sqrt(k*k*(i*i + j*j) + pow(i*i + j*j, 2.));
               nhat[2] = -(i*i + j*j)/sqrt(k*k*(i*i + j*j) + pow(i*i + j*j, 2.)); 
             }
    }

    else if (std::abs(j) > 0) { mhat[0] = 0.; mhat[1] = 0.; mhat[2] = -1.;
                      nhat[0] = -j/sqrt(j*j + i*i);
                      nhat[1] = i/sqrt(j*j + i*i);
                      nhat[2] = 0.; 
                    }

    else if (std::abs(i) > 0) { mhat[0] = 0.; mhat[1] = 1.; mhat[2] = 0.;
                      nhat[0] = 0.; nhat[1] = 0.; nhat[2] = 1.;
                    }

    else if (i==0 && j==0 && k==0) { ; }

    else 
    {
        Error("RandomField::calculate_polarisation_tensors Part of Fourier grid not covered.");
    }

    if (which == 0) { return mhat[l]; }
    else if (which == 1) { return nhat[l]; }
    else { Error("RandomField::calculate_polarisation_tensors Basis vector choice invalid."); }
}

inline void RandomField::init()
{
    BL_PROFILE("RandomField::init_random_field");

    // Set up the problem domain and MF ingredients (Real space)
    IntVect domain_low(0, 0, 0);
    IntVect domain_high(N-1, N-1, N-1);
    Box domain(domain_low, domain_high);
    BoxArray xba(domain);
    DistributionMapping xdm(xba);

    // Make the fft and store the problem domain and MF ingredients (Fourier space)
    FFT::R2C random_field_fft(domain);
    auto const& [kba, kdm] = random_field_fft.getSpectralDataLayout();

    // Set up the arrays to store the in/out data sets
    cMultiFab hs_k(kba, kdm, 2, 0);
    MultiFab hs_x(xba, xdm, 2, 0);
    cMultiFab hij_k(kba, kdm, 6, 0);
    MultiFab hij_x(xba, xdm, 6, 0);

    std::string Filename = "./GRTeclyn-hij-k";
    // Loop to create Fourier-space tensor object
    for (MFIter mfi(hs_k); mfi.isValid(); ++mfi) 
    {
        // Make a pointer to the mode functions at this MF box
        Array4<GpuComplex<Real>> const& hs_ptr = hs_k.array(mfi);
        Array4<GpuComplex<Real>> const& hij_ptr = hij_k.array(mfi);
        const Box& bx = mfi.fabbox();

        // Loop to create mode functions then hij(k)
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            // Find the mode function realisation
            for(int p=0; p<2; p++)
            {
                hs_ptr(i, j, k, p) = calculate_random_field(i, j, k, "position");
            }

            // Find basis tensors and initial tensor realisation
            Vector<Real> eplus(6, 0.);
            Vector<Real> ecross(6, 0.);
            for (int l=0; l<3; l++) for (int p=l; p<3; p++)
            {
                eplus[lut[l][p]] = basis_vector(i, j, k, l, 0)*basis_vector(i, j, k, p, 0) 
                                    - basis_vector(i, j, k, l, 1)*basis_vector(i, j, k, p, 1);
                ecross[lut[l][p]] = basis_vector(i, j, k, l, 0)*basis_vector(i, j, k, p, 1) 
                                    + basis_vector(i, j, k, l, 1)*basis_vector(i, j, k, p, 0);

                hij_ptr(i, j, k, lut[l][p]) = (eplus[lut[l][p]] * hs_ptr(i, j, k, 0)
                                                + ecross[lut[l][p]] * hs_ptr(i, j, k, 1))/std::sqrt(2.);
            }

            // Nyquist node condition
            if ((i==0 || i==N/2) && (j==0 || j==N/2) && (k==0 || k== N/2))
            {
                for(int p=0; p<2; p++) 
                { 
                    GpuComplex<Real> temp(hs_ptr(i, j, k, p).real(), 0.);
                    hs_ptr(i, j, k, p) = temp;
                }

                for (int l=0; l<3; l++) for (int p=l; p<3; p++)
                {
                    GpuComplex<Real> temp(hij_ptr(i, j, k, lut[l][p]).real(), 0.);
                    hij_ptr(i, j, k, lut[l][p]) = temp;
                }
            }

            // Nyquist axis condition
            if (k==0 || k==N/2) 
            {
                if((i>N/2 && j==N/2) || (i==0 && j>N/2) || (i>N/2 && j==0) || (i==N/2 && j>N/2))
                {
                    for(int p=0; p<2; p++) 
                    {
                        GpuComplex<Real> temp(hs_ptr(invert_index(i), invert_index(j), k, p).real(), 
                                                -hs_ptr(invert_index(i), invert_index(j), k, p).imag());
                        hs_ptr(i, j, k, p) = temp;
                    }
                    for (int l=0; l<3; l++) for (int p=l; p<3; p++)
                    {
                        GpuComplex<Real> temp(hij_ptr(invert_index(i), invert_index(j), k, lut[l][p]).real(), 
                                                -hij_ptr(invert_index(i), invert_index(j), k, lut[l][p]).imag());
                        hij_ptr(i, j, k, lut[l][p]) = temp;
                    }
                }
                else if(j > N/2)
                {
                    for(int p=0; p<2; p++) 
                    {
                        GpuComplex<Real> temp(hs_ptr(flip_index(i), invert_index(j), k, p).real(), 
                                                -hs_ptr(flip_index(i), invert_index(j), k, p).imag());
                        hs_ptr(i, j, k, p) = temp;
                    }
                    for (int l=0; l<3; l++) for (int p=l; p<3; p++)
                    {
                        GpuComplex<Real> temp(hij_ptr(flip_index(i), invert_index(j), k, lut[l][p]).real(), 
                                                -hij_ptr(flip_index(i), invert_index(j), k, lut[l][p]).imag());
                        hij_ptr(i, j, k, lut[l][p]) = temp;
                    }
                }
            }

            /*PrintToFile(Filename, 0) << i << "," << j << "," << k;
            for(int s=0; s<6; s++)
            {
                PrintToFile(Filename, 0) << "," << hij_ptr(i, j, k, s) ;
            }
            PrintToFile(Filename, 0) << "\n";*/
        });
    }

    random_field_fft.backward(hs_k, hs_x);
    //random_field_fft.backward(hij_k, hij_x);

    std::string filename = "./GRTeclyn-hij";
    for (MFIter mfi(hs_x); mfi.isValid(); ++mfi) 
    {
        Array4<Real> const& hs_ptr_x = hs_x.array(mfi);
        const Box& bx = mfi.fabbox();

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            PrintToFile(filename, 0) << i << "," << j << "," << k;
            for(int s=0; s<2; s++)
            {
                PrintToFile(filename, 0) << "," << hs_ptr_x(i, j, k, s) ;
            }
            PrintToFile(filename, 0) << "\n";
        });
    }

    //Error("End of first box loop.");
}

#endif /* RANDOMFIELD_IMPL_HPP_*/

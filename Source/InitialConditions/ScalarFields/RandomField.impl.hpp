/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */


#if !defined(RANDOMFIELD_HPP_)
#error "This file should only be included via RandomField.hpp"
#endif

#ifndef RANDOMFIELD_IMPL_HPP_
#define RANDOMFIELD_IMPL_HPP_

inline int RandomField::invert_index(int indx) { return (int)(N/2 - std::abs(N/2 - indx)); }

inline GpuComplex<Real> RandomField::calculate_mode_function(double km, std::string spec_type)
{
    if(km < 1.e-23) { return 0.; }
    
    Real ms_mag = 0.;
    Real ms_arg = 0.;

    Real H0 = sqrt((4.0 * M_PI/3.0/pow(m_params.Mp, 2.))
                * (pow(m_background_params.m * m_background_params.phi0, 2.0) 
                    + pow(m_background_params.Pi0, 2.)));

    double kpr = km/H0;
    if (spec_type == "position")
    {
        ms_mag = sqrt((1.0/km + H0*H0/pow(km, 3.))/2.);
        ms_arg = atan2((cos(kpr) + kpr*sin(kpr)), (kpr*cos(kpr) - sin(kpr)));
    }
    else if (spec_type == "velocity")
    {
        ms_mag = sqrt(km/2.);
        ms_arg = -atan2(cos(kpr), sin(kpr));
    }
    else { Error("RandomField::calculate_power_spectrum Value of spec_type not allowed."); }


    GpuComplex<Real> ps(ms_mag * cos(ms_arg), ms_mag * sin(ms_arg));
    return ps;
}

inline GpuComplex<Real> RandomField::calculate_random_field(int I, int J, int k, std::string spectrum_type)
{
    GpuComplex<Real> value(0., 0.);

    int i = invert_index(I);
    int j = invert_index(J);
    double kmag = std::sqrt(i*i + j*j + k*k) * 2 * M_PI / m_params.L;

    value = calculate_mode_function(kmag, spectrum_type);

    /*if(use_rand)
    {
        double rand_mod = 0.;
        double rand_arg = 0.;
    }*/

    if(m_params.use_window) 
    { 
        double ks = m_params.kstar * 2. * M_PI/m_params.L;
        double Dt = m_params.L/m_params.Delta;
        value *= 0.5 * (1. - tanh(Dt * (kmag - ks))); 
    }

    return value;
}

inline void RandomField::init()
{
    BL_PROFILE("RandomField::init_random_field");
    N = m_params.N_readin;

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
    cMultiFab hij_k(kba, kdm, 6, 0);
    MultiFab hij_x(xba, xdm, 6, 0);

    // Loop to create Fourier-space tensor object
    for (MFIter mfi(hs_k); mfi.isValid(); ++mfi) 
    {
        // Make a pointer to the mode functions at this MF box
        Array4<GpuComplex<Real>> const& hs_ptr = hs_k.array(mfi);
        const Box& bx = mfi.fabbox();

        // Loop to create mode functions then hij(k)
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            hs_ptr(i, j, k, 0) = calculate_random_field(i, j, k, "position");
            hs_ptr(i, j, k, 1) = calculate_random_field(i, j, k, "position");
        });
    }
}

#endif /* RANDOMFIELD_IMPL_HPP_*/

/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */


#if !defined(RANDOMFIELD_HPP_)
#error "This file should only be included via RandomField.hpp"
#endif

#ifndef RANDOMFIELD_IMPL_HPP_
#define RANDOMFIELD_IMPL_HPP_

template <class data_t> void RandomField::init_random_field() const
{
    BL_PROFILE("RandomField::init_random_field");
    N = m_params.N_readin;

    // Set up the problem domain and R2C fft 
    amrex::IntVect domain_low(0, 0, 0);
    amrex::IntVect domain_high(N-1, N-1, N/2);
    amrex::Box domain(domain_low, domain_high);
    //amrex::FFT::R2C random_field_fft(domain);

    // Set up the arrays to store the field data sets
    amrex::BoxArray ba(domain);
    amrex::DistributionMapping dm(ba);
    amrex::MultiFab hs_k(ba, dm, 2, 0);
    amrex::MultiFab hij_k(ba, dm, 6, 0);
    amrex::MultiFab hij_x(ba, dm, 6, 0);
}

#endif /* RANDOMFIELD_IMPL_HPP_*/
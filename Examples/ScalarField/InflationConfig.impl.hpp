/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */


#if !defined(INFLATIONCONFIG_HPP_)
#error "This file should only be included via InflationConfig.hpp"
#endif

#ifndef INFLATIONCONFIG_IMPL_HPP_
#define INFLATIONCONFIG_IMPL_HPP_

// Calculates basis vectors required for polarisation tensors
inline amrex::GpuArray<amrex::Real, 3> 
InflationConfig::calculate_basis_vector(const amrex::IntVect iv, 
                                        const int which_vector)
{
    AMREX_ASSERT(norm > 0);

    // Hermitian symmetry inversion on j and k, with sign on the last two indices.
    // (!!) The FT implemented in AMReX symmetrises across the i index.
    const amrex::Real i = static_cast<amrex::Real>(iv[0]);
    const amrex::Real j = static_cast<amrex::Real>(invert_index_with_sign(iv[1]));
    const amrex::Real k = static_cast<amrex::Real>(invert_index_with_sign(iv[2]));

    amrex::GpuArray<amrex::Real, 3> mhat;
    amrex::GpuArray<amrex::Real, 3> nhat;
    mhat.fill(0.);
    nhat.fill(0.);

    // Skip the 0 mode, as tensors have no average
    if (iv == amrex::IntVect{0, 0, 0}) { ; }

    else if (i > 0.) 
    {
        if (k == 0. && j == 0.) 
        { 
            mhat = amrex::GpuArray<amrex::Real, 3>{0., 1., 0.};
            nhat = amrex::GpuArray<amrex::Real, 3>{0., 0., 1.}; 
        }

        else 
        { 
            amrex::Real norm = sqrt((i*i + j*j) * (i*i + j*j + k*k));
            mhat = amrex::GpuArray<amrex::Real, 3>{j/sqrt(i*i + j*j), -i/sqrt(i*i + j*j), 0.}; 
            nhat = amrex::GpuArray<amrex::Real, 3>{(i*k) / norm, (j*k) / norm, -(i*i + j*j) / norm}; 
        }
    }

    else if (std::abs(j) > 0) 
    { 
        if(k == 0.)
        {
            mhat = amrex::GpuArray<amrex::Real, 3>{0., 0., 1.};
            nhat = amrex::GpuArray<amrex::Real, 3>{1., 0., 0.};
        }

        else
        {
            mhat = amrex::GpuArray<amrex::Real, 3>{-1., 0., 0.};
            nhat = amrex::GpuArray<amrex::Real, 3>{0., -k / sqrt(j*j + k*k), j / sqrt(j*j + k*k)};
        }
    }

    else if (std::abs(k) > 0) 
    { 
        mhat = amrex::GpuArray<amrex::Real, 3>{1., 0., 0.};
        nhat = amrex::GpuArray<amrex::Real, 3>{0., 1., 0.};
    }

    else 
    {
        amrex::Error("RandomField::calculate_polarisation_tensors Part of Fourier grid not covered.");
    }

    if (alpha != 0)
    {
        amrex::Real a = alpha * (M_PI) / 180.;
        amrex::GpuArray<amrex::Real, 3> mp, np;
        for(int l=0; l<3; l++)
        {
            mp[l] = cos(a) * mhat[l] + sin(a) * nhat[l];
            np[l] = -sin(a) * mhat[l] + cos(a) * nhat[l];
        }

        mhat = mp;
        nhat = np;
    }

    if(which_vector == 0) { return mhat; }
    else if(which_vector == 1) { return nhat; }
    else 
    { 
        amrex::Error("RandomField::calculate_basis_vector Incompatable vector type."); 
        return amrex::GpuArray<amrex::Real, 3>{0,0,0}; 
    }
}

// Applies above Nyquist conditions to a given MF
inline void InflationConfig::apply_nyquist_conditions(amrex::cMultiFab &field)
{
    AMREX_ASSERT(N_fine > 0);
    
    int nc = field.nComp();
    for (amrex::MFIter mfi(field); mfi.isValid(); ++mfi) 
    {
        // The geometry for this MPI rank
        const amrex::Box& bx = mfi.fabbox();
        amrex::Array4<amrex::GpuComplex<amrex::Real>> const& field_ptr = field.array(mfi);

        amrex::ParallelFor(bx, [=, this] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            amrex::IntVect iv = {i, j, k};

            if ((i == 0 || i == N_fine/2) && (j == 0 || j == N_fine/2) && (k == 0 || k == N_fine/2))
            {
                for(int comp = 0; comp < nc; comp++)
                {
                    amrex::GpuComplex<amrex::Real> temp(field_ptr(i, j, k, comp).real(), 0.);
                    field_ptr(i, j, k, comp) = temp;
                }
            }

            else if (i==0 || i==N_fine/2) 
            {
                if((k > N_fine/2 && j == N_fine/2) || (k == 0 && j > N_fine/2) ||
                    (k > N_fine/2 && j == 0) || (k == N_fine/2 && j > N_fine/2))
                {
                    for(int comp = 0; comp < nc; comp++) 
                    {
                        amrex::GpuComplex<amrex::Real> temp(field_ptr(i, invert_index(j), invert_index(k), comp).real(), 
                                                -field_ptr(i, invert_index(j), invert_index(k), comp).imag());
                        field_ptr(i, j, k, comp) = temp;
                    }
                }
                
                else if(j > N_fine/2)
                {
                    for(int comp = 0; comp < nc; comp++) 
                    {
                        amrex::GpuComplex<amrex::Real> temp(field_ptr(i, invert_index(j), flip_index(k), comp).real(), 
                                                -field_ptr(i, invert_index(j), flip_index(k), comp).imag());
                        field_ptr(i, j, k, comp) = temp;
                    }
                }
            }
        });

        amrex::Gpu::streamSynchronize();
    }
}

#endif /* INFLATIONCONFIG_IMPL_HPP */
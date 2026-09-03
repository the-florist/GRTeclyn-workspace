/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#if !defined(CONFIG_HPP_)
#error "This file should only be included via Config.hpp"
#endif

#ifndef CONFIG_IMPL_HPP_
#define CONFIG_IMPL_HPP_

// Calculates both basis vectors required for the polarisation tensors
AMREX_GPU_HOST_DEVICE inline Config::BasisVectors
Config::calculate_basis_vectors(const amrex::IntVect iv) const
{
    using Vec = amrex::GpuArray<amrex::Real, 3>;

    // Hermitian symmetry inversion on j and k, with sign on the last two
    // indices. (!!) The FT implemented in AMReX symmetrises across the i index,
    // so i >= 0 always.
    const amrex::Real i = static_cast<amrex::Real>(iv[0]);
    const amrex::Real j =
        static_cast<amrex::Real>(invert_index_with_sign(iv[1]));
    const amrex::Real k =
        static_cast<amrex::Real>(invert_index_with_sign(iv[2]));

    // Default is the zero mode: mhat = nhat = 0, tensors have no average
    Vec mhat{0., 0., 0.};
    Vec nhat{0., 0., 0.};

    if (iv == amrex::IntVect{0, 0, 0})
    { /* zero mode: leave as zero */
    }

    else if (i != 0.)
    {
        if (j == 0. && k == 0.)
        {
            mhat = Vec{0., 1., 0.};
            nhat = Vec{0., 0., 1.};
        }
        else
        {
            const amrex::Real i2j2 = i * i + j * j;
            const amrex::Real ij   = std::sqrt(i2j2);
            const amrex::Real n    = std::sqrt(i2j2 * (i2j2 + k * k));
            mhat                   = Vec{j / ij, -i / ij, 0.};
            nhat                   = Vec{(i * k) / n, (j * k) / n, -i2j2 / n};
        }
    }

    else if (j != 0.) // i == 0
    {
        if (k == 0.)
        {
            mhat = Vec{0., 0., 1.};
            nhat = Vec{1., 0., 0.};
        }
        else
        {
            const amrex::Real jk = std::sqrt(j * j + k * k);
            mhat                 = Vec{-1., 0., 0.};
            nhat                 = Vec{0., -k / jk, j / jk};
        }
    }

    else if (k != 0.) // i == 0, j == 0
    {
        mhat = Vec{1., 0., 0.};
        nhat = Vec{0., 1., 0.};
    }

    else
    {
        // Unreachable: (i, j, k) == (0, 0, 0) only when iv == {0,0,0}, which is
        // handled by the zero-mode branch above.
        AMREX_ASSERT_WITH_MESSAGE(false, "Config::calculate_basis_vectors, "
                                         "Fourier grid point not covered.");
    }

    // Apply the internal rotation in the +/x decomposition basis, if requested
    if (alpha != 0.)
    {
        const amrex::Real a  = alpha * amrex::Math::pi<amrex::Real>() / 180.;
        const amrex::Real ca = std::cos(a);
        const amrex::Real sa = std::sin(a);
        Vec mp, np;
        for (int l = 0; l < 3; l++)
        {
            mp[l] = ca * mhat[l] + sa * nhat[l];
            np[l] = -sa * mhat[l] + ca * nhat[l];
        }
        mhat = mp;
        nhat = np;
    }

    return {mhat, nhat};
}

// Applies above Nyquist conditions to a given MF
inline void Config::apply_nyquist_conditions(amrex::cMultiFab &field)
{
    AMREX_ASSERT(N_fine > 0);

    // Slice to the POD base so the kernel captures config by value, not via
    // the (host) this pointer
    const Config cfg = *this;

    int num_components = field.nComp();
    for (amrex::MFIter mfi(field); mfi.isValid(); ++mfi)
    {
        // The geometry for this MPI rank
        const amrex::Box &bx = mfi.fabbox();
        amrex::Array4<amrex::GpuComplex<amrex::Real>> const &field_ptr =
            field.array(mfi);

        amrex::ParallelFor(
            bx,
            [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
            {
                amrex::IntVect iv = {i, j, k};

                if ((i == 0 || i == cfg.N_fine / 2) &&
                    (j == 0 || j == cfg.N_fine / 2) &&
                    (k == 0 || k == cfg.N_fine / 2))
                {
                    for (int comp = 0; comp < num_components; comp++)
                    {
                        amrex::GpuComplex<amrex::Real> temp(
                            field_ptr(i, j, k, comp).real(), 0.);
                        field_ptr(i, j, k, comp) = temp;
                    }
                }

                else if (i == 0 || i == cfg.N_fine / 2)
                {
                    const int nh = cfg.N_fine / 2;
                    if ((k > nh && j == nh) || (k == 0 && j > nh) ||
                        (k > nh && j == 0) || (k == nh && j > nh))
                    {
                        for (int comp = 0; comp < num_components; comp++)
                        {
                            amrex::GpuComplex<amrex::Real> temp(
                                field_ptr(i, cfg.invert_index(j),
                                          cfg.invert_index(k), comp)
                                    .real(),
                                -field_ptr(i, cfg.invert_index(j),
                                           cfg.invert_index(k), comp)
                                     .imag());
                            field_ptr(i, j, k, comp) = temp;
                        }
                    }

                    else if (j > nh)
                    {
                        for (int comp = 0; comp < num_components; comp++)
                        {
                            amrex::GpuComplex<amrex::Real> temp(
                                field_ptr(i, cfg.invert_index(j),
                                          cfg.flip_index(k), comp)
                                    .real(),
                                -field_ptr(i, cfg.invert_index(j),
                                           cfg.flip_index(k), comp)
                                     .imag());
                            field_ptr(i, j, k, comp) = temp;
                        }
                    }
                }
            });

        amrex::Gpu::streamSynchronize();
    }
}

#endif /* CONFIG_IMPL_HPP_ */
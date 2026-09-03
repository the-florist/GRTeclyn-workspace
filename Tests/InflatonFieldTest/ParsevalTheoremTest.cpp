/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

// AMReX includes
#include <AMReX.H>
#include <AMReX_FFT.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_Math.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Reduce.H>

// Doctest includes
#include "doctest.h"
#include "doctestCLIArgs.hpp"

#include <cmath>

// Test header
#include "ParsevalTheoremTest.hpp"

// For InflatonUtils::tolerance, shared with the rest of the InflatonField
// tests so this test's auto-scaled tolerance stays consistent with theirs.
#include "InflatonUtils.hpp"

namespace
{
// Sum the total power |f_k|^2 over Fourier space, accounting for the
// Hermitian modes not stored explicitly (factor of 2 in the bulk).
amrex::Real calculate_total_power(const amrex::cMultiFab &fk, const int N)
{
    // 1. Set up the parallel reduction operation (Sum)
    amrex::ReduceOps<amrex::ReduceOpSum> reduce_op;
    amrex::ReduceData<amrex::Real> reduce_data(reduce_op);
    using ReduceTuple = typename decltype(reduce_data)::Type;

    // 2. Loop over the grids owned by this MPI rank
    for (amrex::MFIter mfi(fk); mfi.isValid(); ++mfi)
    {
        const amrex::Box &bx = mfi.fabbox();
        auto const &arr      = fk.const_array(mfi);

        // 3. Execute the parallel reduction on the CPU/GPU
        reduce_op.eval(bx, reduce_data,
                       [=] AMREX_GPU_DEVICE(int i, int j, int k) -> ReduceTuple
                       {
                           // Get the real and imaginary parts
                           amrex::Real re = arr(i, j, k).real();
                           amrex::Real im = arr(i, j, k).imag();

                           amrex::Real pow = re * re + im * im;
                           // Multiply by 2 in most of the bulk,
                           // to account for Hermitian modes.
                           if (i != 0 && i != N / 2)
                           {
                               pow *= 2.;
                           }

                           return pow;
                       });
    }

    // 4. Extract the aggregated sum for this specific MPI rank
    ReduceTuple hv                = reduce_data.value();
    amrex::Real total_local_power = amrex::get<0>(hv);

    // 5. Perform an MPI reduction to sum across all compute nodes
    amrex::ParallelDescriptor::ReduceRealSum(total_local_power);

    return total_local_power;
}
} // namespace

void run_parseval_theorem_test()
{
    int amrex_argc    = doctest::cli_args.argc();
    char **amrex_argv = doctest::cli_args.argv();

    // NOLINTNEXTLINE(bugprone-casting-through-void) // Open MPI triggers this
    amrex::Initialize(amrex_argc, amrex_argv);
    {
        constexpr int num_cells = 16;
        const amrex::Real pi    = amrex::Math::pi<amrex::Real>();

        // Config-space (x) and Fourier-space (k) grids, matching the layout
        // used by the InflatonField example's own forward FFTs.
        const amrex::IntVect domain_low(0, 0, 0);
        const amrex::Box x_domain(
            domain_low,
            amrex::IntVect{num_cells - 1, num_cells - 1, num_cells - 1});
        const amrex::BoxArray xba{x_domain};
        const amrex::DistributionMapping xdm{xba};

        const amrex::Box k_domain(
            domain_low,
            amrex::IntVect{num_cells / 2, num_cells - 1, num_cells - 1});
        const amrex::BoxArray kba{k_domain};
        const amrex::DistributionMapping kdm{kba};

        amrex::MultiFab hx{xba, xdm, 1, 0};
        amrex::cMultiFab hk{kba, kdm, 1, 0};
        hk.setVal(0.0);

        // A deterministic, multi-frequency real field with non-trivial power
        // spread across several Fourier modes.
        const auto &hx_arrs = hx.arrays();
        amrex::ParallelFor(
            hx,
            [=] AMREX_GPU_DEVICE(int bx, int i, int j, int k)
            {
                const amrex::Real x_ang = 2. * pi * i / num_cells;
                const amrex::Real y_ang = 2. * pi * j / num_cells;
                const amrex::Real z_ang = 2. * pi * k / num_cells;
                hx_arrs[bx](i, j, k, 0) =
                    std::sin(x_ang) + std::cos(2. * y_ang) +
                    std::sin(3. * z_ang) * std::cos(x_ang);
            });
        amrex::Gpu::streamSynchronize();

        // Forward transform (unnormalised, FFTW-style, as amrex::FFT::R2C
        // returns it) before any physical normalisation is applied.
        amrex::FFT::R2C<amrex::Real> fft(
            x_domain, amrex::FFT::Info().setBatchSize(hk.nComp()));
        fft.forward(hx, hk);

        // Confirm Parseval's theorem holds between config and Fourier space.
        const amrex::Real xsum =
            std::pow(hx.norm2(), 2.) / std::pow(num_cells, 3.);
        const amrex::Real ksum = calculate_total_power(hk, num_cells);

        const int p =
            static_cast<int>(std::round(std::log10((xsum + ksum) / 2.)));
        const amrex::Real tol = InflatonUtils::tolerance * std::pow(10., p + 1);

        INFO("xsum = " << xsum << ", ksum = " << ksum << ", tol = " << tol);
        CHECK(std::abs(xsum - ksum) <= tol);
    }
    amrex::Finalize();
}

/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef TENSORTESTS_HPP_
#define TENSORTESTS_HPP_

#include "InflationUtils.hpp"
#include "Tensor.hpp"

#include <AMReX_MultiFab.H>
#include <AMReX_Reduce.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Print.H>
#include <AMReX_Math.H>

namespace TensorTests
{
    // Test that the input tensor field (config space) is trace free (global)
    inline void Test_is_trace_free(amrex::MultiFab &field)
    {
        if (field.nComp() != 6)
        {
            amrex::Error("RandomField::Test_is_trace_free, "
                         "input field is not a tensor field");
        }

        const auto &arrs = field.arrays();
        amrex::ParallelFor(field,
            [=] AMREX_GPU_DEVICE (int bx, int i, int j, int k)
            {
                amrex::Real sum = 0.;

                for(int l=0; l<3; l++)
                {
                    sum += arrs[bx](i, j, k, InflationUtils::lut[l][l]);
                }

                AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
                    amrex::Math::abs(sum) <= InflationUtils::tolerance,
                    "TensorTests::Test_is_trace_free, trace-free test failed");
            }
        );

        amrex::Gpu::streamSynchronize();
    }

    // Test that the input vectors are orthonormal (local)
    inline void Test_vector_orthonorm(const amrex::IntVect iv,
                                      const amrex::GpuArray<amrex::Real, 3> mhat,
                                      const amrex::GpuArray<amrex::Real, 3> nhat)
    {
        // Confirm basis vectors are orthonormal
        if (iv != amrex::IntVect{0, 0, 0})
        {
            amrex::Real dot1 = 0.;
            amrex::Real dot2 = 0.;
            amrex::Real cross1 = 0.;
            for(int l=0; l<3; l++)
            {
                dot1 += mhat[l] * mhat[l];
                dot2 += nhat[l] * nhat[l];
                cross1 += mhat[l] * nhat[l];
            }

            if (std::abs(dot1 - 1.) > InflationUtils::tolerance 
              || std::abs(dot2 - 1.) > InflationUtils::tolerance 
              || cross1 > InflationUtils::tolerance)
            {
                amrex::Print() << "Location: " << iv << "\n";
                amrex::Print() << "Dot products: " << dot1 << ", " << dot2 << "\n";
                amrex::Print() << "Cross product: " << cross1 << "\n";
                amrex::Print() << "Vectors: \n";
                for(int l=0; l<3; l++)
                {
                    amrex::Print() << l << ", " << mhat[l] << ", " << nhat[l] << "\n";
                }
                amrex::Error("RandomField::Test_vector_orthonorm: "
                             "Basis vectors are not orthonormal here");
            }
        }
    }

    // Test that the input basis tensors, and their rotated counterparts, are orthonormal
    inline void Test_polarisation_tensor_orthonorm(const amrex::IntVect iv,
                                                   const Tensor::Rank2 eplus,
                                                   const Tensor::Rank2 ecross)
    {
        amrex::GpuArray<amrex::Real, 3> conds;
        conds.fill(0.);

        for (int l=0; l<3; l++) for (int p=0; p<3; p++)
        {
            conds[0] += eplus(l, p) * eplus(l, p);
            conds[1] += eplus(l, p) * ecross(l, p);
            conds[2] += ecross(l, p) * ecross(l, p);
        }

        if (iv != amrex::IntVect{0, 0, 0})
        {
            bool plc = (std::abs(conds[0] / 2. - 1.) > InflationUtils::tolerance 
                        || std::abs(conds[2] / 2. - 1.) > InflationUtils::tolerance);
            bool ppc = (std::abs(conds[1]) > InflationUtils::tolerance);
            if (plc || ppc)
            {
                amrex::Print() << "---------\nLocation: " << iv << "\n"
                               << "Dot products: " << conds[0] / 2. << ", "
                               << conds[2] / 2. << "\n"
                               << "Cross product: " << conds[1] << "\n"
                               << "Base tensor components: \n";
                for (int l=0; l<3; l++) for (int p=0; p<3; p++)
                {
                    amrex::Print() << l << ", " << p << ": "
                                   << eplus(l, p) << ", " << ecross(l, p) << "\n";
                }
                amrex::Error("RandomField::Test_polarisation_tensor_orthonorm: "
                             "polarisation tensors are not orthonormal here");
            }
        }
    }

    // Sum the total power |f_k|^2 over Fourier space, accounting for the
    // Hermitian modes not stored explicitly (factor of 2 in the bulk)
    inline amrex::Real calculate_total_power(const amrex::cMultiFab &fk, const int N)
    {
        // 1. Set up the parallel reduction operation (Sum)
        amrex::ReduceOps<amrex::ReduceOpSum> reduce_op;
        amrex::ReduceData<amrex::Real> reduce_data(reduce_op);
        using ReduceTuple = typename decltype(reduce_data)::Type;

        // 2. Loop over the grids owned by this MPI rank
        for (amrex::MFIter mfi(fk); mfi.isValid(); ++mfi)
        {
            const amrex::Box &bx = mfi.fabbox();
            auto const &arr = fk.const_array(mfi);

            // 3. Execute the parallel reduction on the CPU/GPU
            reduce_op.eval(bx, reduce_data,
                [=] AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
                {
                    // Get the real and imaginary parts
                    amrex::Real re = arr(i, j, k).real();
                    amrex::Real im = arr(i, j, k).imag();

                    amrex::Real pow = re * re + im * im;
                    // Multiply by 2 in most of the bulk,
                    // to account for Hermitian modes.
                    if (i != 0 && i != N/2) { pow *= 2.; }

                    return pow;
                });
        }

        // 4. Extract the aggregated sum for this specific MPI rank
        ReduceTuple hv = reduce_data.value();
        amrex::Real total_local_power = amrex::get<0>(hv);

        // 5. Perform an MPI reduction to sum across all compute nodes
        amrex::ParallelDescriptor::ReduceRealSum(total_local_power);

        return total_local_power;
    }

    // Confirm Parseval's theorem holds between a config-space field hx and its
    // Fourier-space counterpart hk (checked before physical normalisation)
    inline void Test_Parsevals_thm(const amrex::MultiFab &hx,
                                   const amrex::cMultiFab &hk, const int N)
    {
        amrex::Real xsum = std::pow(hx.norm2(), 2.);
        xsum /= std::pow(N, 3.);

        amrex::Real ksum = calculate_total_power(hk, N);

        int p = std::round(std::log10((ksum + ksum) / 2.));
        amrex::Real tol = InflationUtils::tolerance * std::pow(10., p+1);

        if (std::abs(xsum - ksum) > tol)
        {
            amrex::Print() << "Order of magnitude of sum: " << p << "\n";
            amrex::Print() << "Tolerance: " << tol << "\n";
            amrex::Print() << "Stdev (x): " << xsum << "\n";
            amrex::Print() << "Integrated power (k): " << ksum << "\n";
            amrex::Print() << "Ratio: " << ksum / xsum << "\n";
            amrex::Print() << "Difference: " << std::abs(ksum - xsum) << "\n";
            amrex::Error("TensorTests::Test_Parsevals_thm, "
                         "Parseval's theorem fails here.");
        }
    }
}

#endif /* TENSORTESTS_HPP_ */
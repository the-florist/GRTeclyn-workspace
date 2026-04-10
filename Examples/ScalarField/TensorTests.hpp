/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef TENSORTESTS_HPP_
#define TENSORTESTS_HPP_

#include "InflationUtils.hpp"

namespace TensorTests 
{
    // Test that the input tensor field (config space) is trace free (global)
    inline void Test_is_trace_free(amrex::MultiFab &field)
    {
        if (field.nComp() != 6)
        {
            amrex::Error("RandomField::Test_is_trace_free, input field is not a tensor field");
        }

        const auto &arrs = field.arrays();
        amrex::ParallelFor(field,
            [=] AMREX_GPU_DEVICE (int bx, int i, int j, int k)
            {
                amrex::IntVect iv{i, j, k};
                amrex::Real sum = 0.;

                for(int l=0; l<3; l++)
                {
                    sum += arrs[bx](i, j, k, InflationUtils::lut[l][l]);
                }

                if(std::abs(sum) > InflationUtils::tolerance)
                {
                    amrex::Print() << iv << ": " << sum << "\n";
                    amrex::Error("RandomField::Test_is_trace_free Trace-free test failed here.");
                }
            }
        );
    }

    // Test that the input vectors are orthonormal (local)
    inline void Test_vector_orthonorm(const amrex::IntVect iv, const amrex::Vector<amrex::Real> mhat, 
                                            const amrex::Vector<amrex::Real> nhat)
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

            if(std::abs(dot1 - 1.) > InflationUtils::tolerance 
              || std::abs(dot2 - 1.) > InflationUtils::tolerance 
              || cross1 > InflationUtils::tolerance)
            {
                amrex::Print() << "Location: " << iv << "\n";
                amrex::Print() << "Dot products: " << dot1 << ", " << dot2 << "\n";
                amrex::Print() << "Cross product: " << cross1 << "\n";
                amrex::Print() << "amrex::Vectors: \n";
                for(int l=0; l<3; l++)
                {
                    amrex::Print() << l << ", " << mhat[l] << ", " << nhat[l] << "\n";
                }
                amrex::Error("RandomField::Test_vector_orthonorm: Basis vectors are not orthonormal here");
            }
        }
    }

    // Test that the input basis tensors, and their rotated counterparts, are orthonormal
    inline void Test_polarisation_tensor_orthonorm(const amrex::IntVect iv, const Tensor<2, amrex::Real> eplus,
                                                                    const Tensor<2, amrex::Real> ecross)
    {
        amrex::Vector<amrex::Real> conds(3, 0.);

        for (int l=0; l<3; l++) for (int p=0; p<3; p++)
        {
            conds[0] += eplus[l][p] * eplus[l][p];
            conds[1] += eplus[l][p] * ecross[l][p];
            conds[2] += ecross[l][p] * ecross[l][p];
        }

        if(iv != amrex::IntVect{0, 0, 0})
        {
            bool plc = (std::abs(conds[0] / 2. - 1.) > InflationUtils::tolerance 
                        || std::abs(conds[2] / 2. - 1.) > InflationUtils::tolerance);
            bool ppc = (std::abs(conds[1]) > InflationUtils::tolerance);
            if (plc || ppc)
            {
                amrex::Print() << "---------\nLocation: " << iv << "\n";
                amrex::Print() << "Dot products: " << conds[0] / 2. << ", " << conds[2] / 2. << "\n";
                amrex::Print() << "Cross product: " << conds[1] << "\n";
                amrex::Print() << "Base tensor components: \n";
                for (int l=0; l<3; l++) for (int p=0; p<3; p++)
                {
                    amrex::Print() << l << ", " << p << ": ";
                    amrex::Print() << eplus[l][p] << ", " << ecross[l][p] << "\n";
                }
                amrex::Error("RandomField::Test_polarisation_tensor_orthonorm: polarisation tensors are not orthonormal here");
            }
        }
    }
}

#endif /* TENSORTESTS_HPP_ */
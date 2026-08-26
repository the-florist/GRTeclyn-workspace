/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef INFLATIONUTILS_HPP_
#define INFLATIONUTILS_HPP_

#include <cstdint>

#include <AMReX_GpuQualifiers.H>
#include <AMReX_REAL.H>

namespace InflationUtils
{
    // Look-up table
    // Used to construct polarisation basis tensors
    const Tensor<2, int> lut{0, 1, 2, 1, 3, 4, 2, 4, 5};
    constexpr amrex::Real tolerance = 1.e-12;

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    std::uint64_t splitmix64(std::uint64_t x)
    {
        x += 0x9E3779B97F4A7C15ULL;
        x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ULL;
        x = (x ^ (x >> 27)) * 0x94D049BB133111EBULL;
        return x ^ (x >> 31);
    }

    // Map 64 random bits to the OPEN interval (0, 1), so that the Rayleigh
    // draw sqrt(-2 ln u) stays finite (u is never exactly 0 or 1).
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real to_unit_open(const std::uint64_t bits)
    {
        return (static_cast<amrex::Real>(bits >> 11) + 0.5) * (1.0 / 9007199254740992.0);
    }
}

#endif /* INFLATIONUTILS_HPP_ */
/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef INFLATONPARAMETERS_HPP_
#define INFLATONPARAMETERS_HPP_

#include "GRParmParse.hpp"

#include <AMReX_Math.H>
#include <AMReX_Vector.H>

#include "Potential.hpp"

struct InflatonParameters
{
    // Parameters, populated once on the host by fill_params()
    amrex::Real Mp{1.};    //!< Planck mass
    amrex::Real L{0.};     //!< Grid length [Mp^(-1)]
    amrex::Real Delta{1.}; //!< Window function width
    amrex::Real alpha{0.}; //!< Polarisation basis angle

    int N{0};                 //!< Grid resolution
    int N_fine{0};            //!< Finest resolution (convergence tests only)
    int N_coarse{0};          //!< Coarsest resolution (sets cut-off mode)
    int random_seed{3539263}; //!< Seed for random number generator

    int scalar_init{0}; //!< Add scalar perturbations (1) or not (0)
    int tensor_init{0}; //!< Add tensor perturbations (1) or not (0)
    int use_rand{1};    //!< Make perturbations Gaussian random fields
    int use_window{0};  //!< Apply window function to initial spectrum

    amrex::Real phi0{0.};         //!< Background scalar field value
    amrex::Real Pi0{0.};          //!< Background Pi value
    amrex::Real V_background{0.}; //!< Background potential value
    amrex::Real dV_backgrond{0.}; //!< Background potential phi derivative
    amrex::Real H0{0.};           //!< Initial Hubble parameter
    amrex::Real epsilon_1{0.};    //!< First slow-roll parameter
    amrex::Real epsilon_2{0.};    //!< Second slow-roll parameter
    amrex::Real init_a{1.};       //!< Initial scale factor (almost always 1)

    // Read the scalar parameters once from the global GRParmParse table
    void fill_params()
    {
        // Read grid parameters
        GRParmParse pp;
        amrex::Vector<int> ncell;
        pp.getarr("amr.n_cell", ncell);
        N      = ncell[0];
        N_fine = N;
        pp.query("N_fine", N_fine);
        N_coarse = N;
        pp.query("N_coarse", N_coarse);
        amrex::Vector<amrex::Real> prob_extent;
        pp.getarr("geometry.prob_extent", prob_extent);
        L = prob_extent[0];

        // Read and set physical units
        amrex::Real G_Newton = 1.;
        pp.query("scalar_field.G_Newton", G_Newton);
        Mp = 1. / std::sqrt(G_Newton);

        // Read initialisation parameters
        GRParmParse init_pp("init");
        init_pp.query("scalar_init", scalar_init);
        init_pp.query("tensor_init", tensor_init);
        init_pp.query("use_rand", use_rand);
        init_pp.query("use_window", use_window);
        init_pp.query("alpha", alpha);
        init_pp.query("Delta", Delta);
        init_pp.query("random_seed", random_seed);
        init_pp.get("background_phi", phi0);
        init_pp.query("background_dphi", Pi0);

        // Set background spacetime according to Friedmann equations
        Potential potential;
        potential.compute_background_potential(V_background, dV_backgrond,
                                               phi0);
        H0 = calculate_H0(G_Newton, Pi0, V_background);

        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
            scalar_init == 0 || Pi0 != 0.,
            "Config::fill_params, init.background_dphi must be "
            "non-zero when init.scalar_init is enabled");
        if (Pi0 != 0.)
        {
            epsilon_1 = std::pow(Pi0 / H0, 2.) / 2. / std::pow(Mp, 2.);
            epsilon_2 =
                6. + dV / Pi0 / H0 - std::pow(Pi0 / H0, 2.) / std::pow(Mp, 2.);
        }
    }

    void check_params(const amrex::Vector<int> &ncell,
                      const amrex::Vector<amrex::Real> &prob_extent) const
    {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
            Mp != 0., "Config::check_params, Mp must be non-zero");
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
            N_fine >= N, "Config::check_params, N_fine must be >= N");
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
            N_coarse <= N, "Config::check_params, N_coarse must be <= N");
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
            L != 0., "Config::check_params, L must be non-zero");
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
            Delta != 0., "Config::check_params, Delta must be non-zero");

        for (int d = 1; d < static_cast<int>(ncell.size()); d++)
        {
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
                ncell[d] == ncell[0],
                "Config::check_params, grid must be cubic "
                "(amr.n_cell must be equal in all directions)");
        }
        for (int d = 1; d < static_cast<int>(prob_extent.size()); d++)
        {
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
                prob_extent[d] == prob_extent[0],
                "Config::check_params, domain must be cubic "
                "(geometry.prob_extent must be equal in all directions)");
        }
    }

    // Helper function that stores our form of the Friedmann equations
    AMREX_GPU_HOST_DEVICE static amrex::Real
    calculate_H0(amrex::Real G, amrex::Real Pi, amrex::Real V)
    {
        return sqrt((8. * amrex::Math::pi<amrex::Real>() * G / 3.) *
                    (0.5 * pow(Pi, 2.) + V));
    }
};

#endif /* INFLATONPARAMETERS_HPP_ */
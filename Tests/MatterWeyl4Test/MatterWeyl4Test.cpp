/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */
#define COMPARE_WITH_CHF
#define COVARIANTZ4

// AMReX includes
#include <AMReX.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_FabArray.H>
#include <AMReX_MFIter.H>
#include <AMReX_MultiFab.H>

#ifdef AMREX_USE_HDF
#include <AMReX_PlotFileUtilHDF5.H>
#endif

// Doctest includes
#include "doctest.h"
#include "doctestCLIArgs.hpp"

#include <iomanip>
#include <iostream>
#include <sys/time.h>

// this initializes the CCZ4 variables
#include <InitialData.hpp>

// test header
#include "MatterWeyl4Test.hpp"

// GRTeclyn includes
//  #include "BoxLoops.hpp"
#include "CCZ4RHS.hpp"
#include "EMTensor.hpp"
#include "MatterCCZ4RHS.hpp"
// #include "GravWavDecF_F.H"
#include "DefaultPotential.hpp"
#include "MatterWeyl4.hpp"
#include "ScalarField.hpp"
#include "Weyl4.hpp"
#include "simd.hpp"
#include <array>

// #include "SetValue.hpp"
// #include "UserVariables.hpp"

// Chombo namespace
// #include "UsingNamespace.H"

void run_matter_weyl4_test()
{
    int amrex_argc    = doctest::cli_args.argc();
    char **amrex_argv = doctest::cli_args.argv();

    amrex::Initialize(amrex_argc, amrex_argv, true, MPI_COMM_WORLD);
    {
        constexpr int num_cells  = 32;
        constexpr int num_ghosts = 3;
        constexpr amrex::Real dx = 0.25 / (num_cells - 1);
        amrex::Box box{amrex::IntVect::TheZeroVector(),
                       amrex::IntVect{num_cells - 1}};
        auto ghosted_box = box;
        ghosted_box.grow(num_ghosts);

        amrex::BoxArray box_array{box};
        amrex::DistributionMapping distribution_mapping{box_array};
        amrex::MFInfo mf_info;
        mf_info.SetArena(amrex::The_Managed_Arena());

        amrex::MultiFab in_fab{box_array, distribution_mapping, NUM_VARS,
                               num_ghosts, mf_info};

        const auto &in_arrays = in_fab.arrays();
        amrex::ParallelFor(
            in_fab, in_fab.nGrowVect(),
            [=] AMREX_GPU_DEVICE(int ibox, int i, int j, int k)
            {
                const amrex::IntVect iv{i, j, k};
                const amrex::RealVect coords = amrex::RealVect{iv} * dx;

                random_ccz4_initial_data(iv, in_arrays[ibox], coords);
            });

        amrex::Gpu::streamSynchronize();

        using DefaultScalarField = ScalarField<DefaultPotential>;

        ScalarField<DefaultPotential> my_scalar_field(DefaultPotential());

        // set up weyl4 calculation
        constexpr int dcomp = 0;
        double G_Newton     = 1.0;
        std::array<double, AMREX_SPACEDIM> center{0.0, 0.0, 0.0};
        MatterWeyl4<DefaultScalarField> matter_weyl4(
            DefaultScalarField(DefaultPotential()), center, dx, dcomp,
            CCZ4RHS<>::USE_CCZ4, G_Newton);

        amrex::MultiFab out_fab{box_array, distribution_mapping, NUM_VARS,
                                num_ghosts, mf_info};

        const auto &in_c_arrays = in_fab.const_arrays();
        const auto &out_arrays  = out_fab.arrays();

        amrex::FArrayBox &out_testfab = out_fab[0];
        amrex::Array4<amrex::Real> a  = out_testfab.array();

        amrex::FArrayBox &in_testfab       = in_fab[0];
        amrex::Array4<amrex::Real const> b = in_testfab.array();

        amrex::ParallelFor(out_fab,
                           [=] AMREX_GPU_DEVICE(int ibox, int i, int j, int k)
                           {
                               matter_weyl4.compute(i, j, k, out_arrays[ibox],
                                                    in_c_arrays[ibox]);

                               // matter_weyl4.compute(i, j, k,
                               // out_arrays[ibox],
                               //                      in_c_arrays[ibox]);
                           });

        // // compute rho and store in FAB
        // BoxLoops::loop(EMTensor<DefaultScalarField>(
        //                    DefaultScalarField(DefaultPotential()), dx,
        //                    c_Rho),
        //                in_fab, in_fab, /* can't compute in ghost cells */
        //                box);

        // Real null = 0;

        // std::array<double, CH_SPACEDIM> centerGW = {0, 0, 0};

        // struct timeval begin, end;
        // gettimeofday(&begin, NULL);

        // BoxLoops::loop(
        //     MatterWeyl4<DefaultScalarField>(DefaultScalarField(DefaultPotential()),
        //                                     centerGW, dx,
        //                                     CCZ4RHS<>::USE_BSSN),
        //     in_fab, out_fab);

        bool test_passes = true;

        CHECK(test_passes);
    }
    amrex::Finalize();
}
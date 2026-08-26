/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef SIMULATIONPARAMETERS_HPP_
#define SIMULATIONPARAMETERS_HPP_

// General includes
#include "GRParmParse.hpp"
#include "SimulationParametersBase.hpp"

// Problem specific includes:
#include "InitialBackgroundData.hpp"
#include "InitialScalarData.hpp"
#include "Potential.hpp"
#include "InflationConfig.hpp"

class SimulationParameters : public SimulationParametersBase
{
  public:
    SimulationParameters(GRParmParse &pp) : SimulationParametersBase(pp)
    {
        // read the problem specific params
        read_params(pp);
        check_params();
    }

    void read_params(GRParmParse &pp)
    {
        // GR and default parameters
	    initial_params.center =
            center; // already read in SimulationParametersBase
         pp.load("G_Newton", G_Newton,
                 0.0); // for now the example neglects backreaction
        inflt_params.Mp = 1./std::sqrt(G_Newton);
        background_params.G_Newton = G_Newton;

        //<! Random Field init parameters
        // Flags
        pp.load("randominit.read_from_STOIIC", inflt_params.read_from_stoiic, 0);
        pp.load("randominit.tensor_init", inflt_params.tensor_init, 0);
        pp.load("randominit.scalar_init", inflt_params.scalar_init, 0);
        pp.load("randominit.use_rand", inflt_params.use_rand, 1);
        pp.load("randominit.use_window", inflt_params.use_window, 0);

        // Grid parameters
        pp.load("L", inflt_params.L, 1.);
        pp.load("N", inflt_params.N, 32);
        pp.load("N_fine", inflt_params.N_fine, inflt_params.N);
        pp.load("N_coarse", inflt_params.N_coarse, inflt_params.N);

        // Field construction parameters
        pp.load("randominit.A", inflt_params.A, 1.);
        pp.load("randominit.random_seed", inflt_params.random_seed, 3539263);
        pp.load("randominit.alpha", inflt_params.alpha, 0.);
        pp.load("randominit.kstar", inflt_params.kstar, 0.);
        pp.load("randominit.Delta", inflt_params.Delta, 1.);

        // Extraction parameters
        pp.load("extraction.calc_binned_power_spectrum", 
                inflt_params.calc_binned_power_spectrum, 0);
        pp.load("extraction.spec_interval", inflt_params.plot_int, 100);
        pp.load("extraction.calc_higher_order_statistics", 
                inflt_params.calc_higher_order_statistics, 0); 
        pp.load("extraction.num_moments", inflt_params.num_orders, 0);
        pp.getarr("extraction.moments_to_print", inflt_params.orders, 0, 
                                                inflt_params.num_orders);	

        //<! If reading in from STOIIC, read in initial spectrum arrays
        if(inflt_params.read_from_stoiic)
        {
            int num_modes;
            pp.load("randominit.n_k", num_modes, 0);
            pp.getarr("randominit.init_k", inflt_params.init_k, 0, num_modes);

            amrex::Print() << "Begin read in of scalars...\n";

            if(inflt_params.scalar_init)
            {
                inflt_params.scalar_ps = amrex::Vector<amrex::Vector<amrex::Real>>(8, 
                                            amrex::Vector<amrex::Real>(num_modes, 0.));
                pp.getarr("randominit.re_phi_k", inflt_params.scalar_ps[0], 0, num_modes);
                pp.getarr("randominit.im_phi_k", inflt_params.scalar_ps[1], 0, num_modes);
                pp.getarr("randominit.re_Pi_k", inflt_params.scalar_ps[2], 0, num_modes);
                pp.getarr("randominit.im_Pi_k", inflt_params.scalar_ps[3], 0, num_modes);
                pp.getarr("randominit.re_X_k", inflt_params.scalar_ps[4], 0, num_modes);
                pp.getarr("randominit.im_X_k", inflt_params.scalar_ps[5], 0, num_modes);
                pp.getarr("randominit.re_K_k", inflt_params.scalar_ps[6], 0, num_modes);
                pp.getarr("randominit.im_K_k", inflt_params.scalar_ps[7], 0, num_modes);
            }

            if(inflt_params.tensor_init == 1)
            {
                amrex::Print() << "Begin read in of tensors...\n";

                inflt_params.tensor_ps = amrex::Vector<amrex::Vector<amrex::Real>>(4, 
                                            amrex::Vector<amrex::Real>(num_modes, 0.));
                pp.getarr("randominit.re_h_k", inflt_params.tensor_ps[0], 0, num_modes);
                pp.getarr("randominit.im_h_k", inflt_params.tensor_ps[1], 0, num_modes);
                pp.getarr("randominit.re_dh_k", inflt_params.tensor_ps[2], 0, num_modes);
                pp.getarr("randominit.im_dh_k", inflt_params.tensor_ps[3], 0, num_modes);
            }
        }

        //<! Potential parameters
        pp.load("potential.type", potential_params.type, 0);
        pp.load("potential.param_1", potential_params.param1, 0.1);
        pp.load("potential.param_2", potential_params.param2, 0.);
        pp.load("potential.param_3", potential_params.param3, 0.);
        pp.load("potential.param_4", potential_params.param4, 0.);
        pp.load("potential.param_5", potential_params.param5, 0.);

        //<! Background parameters
        pp.load("init.background_phi", background_params.phi0, 0.0);
        pp.load("init.background_dphi", background_params.Pi0, 0.0);
    }

    void check_params()
    {
        warn_parameter("kstar", inflt_params.kstar,
                       inflt_params.kstar >= 0,
                       "cut-off frequency index must be positive");

        check_parameter("Delta", inflt_params.Delta,
                       (!inflt_params.calc_binned_power_spectrum
                        || inflt_params.Delta > 0),
                       "cut-off width must be positive and non-zero");

        check_parameter("orders", inflt_params.calc_higher_order_statistics,
                       (!inflt_params.calc_higher_order_statistics 
                        || !inflt_params.orders.empty()),
                       "moment orders must be provided");

        check_parameter("N_fine", inflt_params.N_fine, 
                        N_fine >= N,
                        "finest resolution must be larger than N");

        warn_parameter("N_coarse", inflt_params.N_coarse, 
                       N_coarse <= N,
                       "coarsest resolution should be smaller than N");
    }

    // Initial data for matter and potential and BH
    double G_Newton;
    Potential::params_t potential_params;
    InitialBackgroundData::params_t background_params;
    InitialScalarData::params_t initial_params;
    InflationConfig inflt_params;
};

#endif /* SIMULATIONPARAMETERS_HPP_ */

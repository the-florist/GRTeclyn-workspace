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
	    initial_params.center =
            center; // already read in SimulationParametersBase
         pp.load("G_Newton", G_Newton,
                 0.0); // for now the example neglects backreaction
        pp.load("scalar_amplitude", initial_params.amplitude, 0.1);
        pp.load("scalar_width", initial_params.width, 1.0);
	pp.load("scalar_mass", potential_params.scalar_mass, 0.1);	

        // Initial scalar field data
        pp.load("G_Newton", G_Newton,
                0.0); // for now the example neglects backreaction
	initial_params.center =
            center; // already read in SimulationParametersBase
	pp.load("scalar_amplitude", background_params.phi0, 0.0);
	pp.load("scalar_velocity", background_params.Pi0, 0.0);
        pp.load("scalar_mass", background_params.m, 0.0);
    }

    void check_params()
    {
        warn_parameter("scalar_mass", background_params.m,
                       background_params.m <
                           0.2 / coarsest_dx / dt_multiplier,
                       "oscillations of scalar field do not appear to be "
                       "resolved on coarsest level");
    }

    // Initial data for matter and potential and BH
    double G_Newton;
    Potential::params_t potential_params;
    InitialBackgroundData::params_t background_params;
    InitialScalarData::params_t initial_params;
};

#endif /* SIMULATIONPARAMETERS_HPP_ */

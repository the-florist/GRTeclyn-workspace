/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef SIMULATIONPARAMETERS_HPP_
#define SIMULATIONPARAMETERS_HPP_

#include "Config.hpp"
#include "Potential.hpp"

class SimulationParameters
{
  public:
    SimulationParameters() = delete;

    static void check_params()
    {
        Config::params_t::check_params();
        Potential::params_t::check_params();
    }
};

#endif /* SIMULATIONPARAMETERS_HPP_ */

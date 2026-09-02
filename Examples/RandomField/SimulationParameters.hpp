/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef SIMULATIONPARAMETERS_HPP_
#define SIMULATIONPARAMETERS_HPP_

#include "BaseParameterChecker.hpp"
#include "CCZ4RHS.hpp"
#include "FixedGridsTagger.hpp"
#include "IntegratedMovingPunctureGauge.hpp"
#include "LineExtractionParameters.hpp"
#include "Potential.hpp"
#include "ScalarField.hpp"

class SimulationParameters
{
  public:
    SimulationParameters() = delete;

    static void check_params()
    {
        ;
    }
};

#endif /* SIMULATIONPARAMETERS_HPP_ */

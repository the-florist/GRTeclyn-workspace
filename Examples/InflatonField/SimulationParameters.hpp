/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef SIMULATIONPARAMETERS_HPP_
#define SIMULATIONPARAMETERS_HPP_

#include "GRParmParse.hpp"
#include "InflatonUtils.hpp"
#include "Potential.hpp"

class SimulationParameters
{
  public:
    SimulationParameters() = delete;

    static void check_params()
    {
        GRParmParse pp;
        amrex::Vector<int> ncell;
        pp.getarr("amr.n_cell", ncell);
        amrex::Vector<amrex::Real> prob_extent;
        pp.getarr("geometry.prob_extent", prob_extent);

        const InflatonUtils utils; // fill_params() runs in the constructor
        utils.m_params.check_params(ncell, prob_extent);

        Potential::params_t::check_params();
    }
};

#endif /* SIMULATIONPARAMETERS_HPP_ */

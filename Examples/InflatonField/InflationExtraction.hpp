/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef INFLATIONEXTRACTION_HPP_
#define INFLATIONEXTRACTION_HPP_

#include "DerivedVariables.hpp"
#include "FilesystemTools.hpp"
#include "GRParmParse.hpp"
#include "InflatonUtils.hpp"
#include "SmallDataIO.hpp"
#include "StateVariables.hpp"

#include <AMReX_MultiFab.H>
#include <AMReX_Vector.H>

#include <string>

//! Diagnostics written to data files during the evolution.
/*!
    Two products are written, both derived from the comoving curvature
    perturbation R and the tensor polarisation fields h+ and hx:

      - binned isotropic power spectra of each field, one file per output
        time, written every `extraction.spec_interval` steps;
      - the first four statistical moments of each field, appended to a
        single time series.

    The configuration- and Fourier-space fields themselves come from
    DerivedVariables::extract_hs_and_R, which owns the transform pipeline;
    this class only bins and reduces them. The plotfile output of the same
    fields is handled separately by DerivedVariables as an AMReX derive.
*/
class InflationExtraction
{
  public:
    //! Field names, in the component order used by extract(). Matches the
    //! ordering of the MultiFab assembled there: R first, then the
    //! polarisations.
    static inline const amrex::Vector<std::string> var_names{"R", "hplus",
                                                             "hcross"};

    struct params_t
    {
        //! Write the volume-averaged background quantities every step
        int calc_background_means{1};

        //! Write binned isotropic power spectra
        int calc_binned_power_spectrum{0};

        //! Number of steps between power spectrum outputs
        int spec_interval{100};

        //! Write the statistical moments time series
        int calc_higher_order_statistics{0};

        //! Which moments to write, any subset of {1, 2, 3, 4} in ascending
        //! order: mean, standard deviation, skewness, kurtosis
        amrex::Vector<int> moment_orders;

        //! Directory that the data files are written into
        std::string data_path{"data/"};

        void fill_params()
        {
            GRParmParse extraction_pp("extraction");
            extraction_pp.query("calc_background_means", calc_background_means);
            extraction_pp.query("calc_binned_power_spectrum",
                                calc_binned_power_spectrum);
            extraction_pp.query("spec_interval", spec_interval);
            extraction_pp.query("calc_higher_order_statistics",
                                calc_higher_order_statistics);
            extraction_pp.query("path", data_path);

            if (extraction_pp.contains("moments_to_print"))
            {
                extraction_pp.getarr("moments_to_print", moment_orders);
            }
        }

        static void check_params()
        {
            GRParmParse extraction_pp("extraction");

            int calc_higher_order_statistics{0};
            extraction_pp.queryAdd("calc_higher_order_statistics",
                                   calc_higher_order_statistics);

            amrex::Vector<int> moment_orders;
            if (extraction_pp.contains("moments_to_print"))
            {
                extraction_pp.getarr("moments_to_print", moment_orders);
            }

            if (calc_higher_order_statistics != 0 && moment_orders.empty())
            {
                extraction_pp.error("moments_to_print",
                                    "must list the moments to write when "
                                    "calc_higher_order_statistics is enabled");
            }

            int previous_order = 0;
            for (const int order : moment_orders)
            {
                if (order < 1 || order > 4)
                {
                    extraction_pp.error("moments_to_print",
                                        "only moments 1 to 4 (mean, stdev, "
                                        "skewness, kurtosis) are implemented");
                }
                if (order <= previous_order)
                {
                    extraction_pp.error(
                        "moments_to_print",
                        "must be given in strictly ascending order");
                }
                previous_order = order;
            }

            int spec_interval{100};
            extraction_pp.queryAdd("spec_interval", spec_interval);
            if (spec_interval <= 0)
            {
                extraction_pp.error("spec_interval", "must be > 0");
            }
        }
    };

    // NOLINTBEGIN(bugprone-easily-swappable-parameters)
    InflationExtraction(params_t a_params, const amrex::Real a_dt,
                        const amrex::Real a_time,
                        const amrex::Real a_restart_time,
                        const bool a_first_step)
        : m_params(std::move(a_params)), m_dt(a_dt), m_time(a_time),
          m_restart_time(a_restart_time), m_first_step(a_first_step)
    {
    }
    // NOLINTEND(bugprone-easily-swappable-parameters)

    //! Writes the requested spectra and statistics for this time step.
    void extract(const amrex::MultiFab &state);

    //! Writes the requested moments of an arbitrary MultiFab as one time-data
    //! line, returning the per-component standard deviations. Public because
    //! the level class also uses it for the constraint variables.
    amrex::Vector<amrex::Real>
    print_moment(const amrex::MultiFab &field,
                 const amrex::Vector<std::string> &names,
                 const amrex::Vector<int> &moment_orders, SmallDataIO &file);

  private:
    params_t m_params;
    amrex::Real m_dt;
    amrex::Real m_time;
    amrex::Real m_restart_time;
    bool m_first_step;

    //! Grid geometry, needed to map a Fourier index onto |k|
    InflatonUtils m_utils;

    [[nodiscard]] const InflatonParameters &params() const
    {
        return m_utils.m_params;
    }

    //! Creates, and returns the path to, a subdirectory of the data directory
    [[nodiscard]] std::string
    make_subdirectory(const std::string &subdirectory) const;

    //! Volume averages of the evolved background quantities, appended to
    //! means-file.dat as one line per time step
    void extract_background_means(const amrex::MultiFab &state);

    //! Bins one component of a Fourier-space field onto an isotropic k axis
    //! and writes the result
    void print_power_spectrum(const amrex::cMultiFab &field_array,
                              SmallDataIO &power_spec_file, int component);

    //! Central moment of the requested order about \p mean
    [[nodiscard]] amrex::Real
    calculate_field_moment_x(const amrex::MultiFab &field, amrex::Real mean,
                             int moment, int component) const;

    //! Places one statistic into the flat header/data layout used by the
    //! statistics file
    static void assign_statistics_data(
        amrex::Vector<std::string> &header_storage, const std::string &name,
        amrex::Vector<amrex::Real> &data_storage, amrex::Real value,
        int component, int num_comps, amrex::Vector<int>::const_iterator itr,
        amrex::Vector<int>::const_iterator start, bool is_first_step);
};

#include "InflationExtraction.impl.hpp"

#endif /* INFLATIONEXTRACTION_HPP_ */

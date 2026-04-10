/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef INFLATIONEXTRACTION_HPP_
#define INFLATIONEXTRACTION_HPP_

#include "InflationConfig.hpp"
#include "VarsTools.hpp"
#include "FilesystemTools.hpp" 
#include "SmallDataIO.hpp"

#include <AMReX_Print.H>
#include <AMReX_Vector.H>
#include <AMReX_FFT.H>

class InflationExtraction
{
    public:
        // Names of diagnostic variables
        static inline const amrex::GpuArray<std::string, 3> var_names = {"R", "hplus", "hcross"};

        // Constructor used in extraction of diagnostics
        InflationExtraction(InflationConfig a_params)
                : m_params(a_params), norm(std::pow(a_params.L, -3.))   
        {}

        // Main routines
        void derive(const amrex::MultiFab &source, amrex::MultiFab &out, const int dcomp);
        void extract(const amrex::MultiFab &state);

        // Set-up function when printing to data files
        void set_print_params(const std::string data_path, const amrex::Real dt,  
                     const amrex::Real cur_time, const int restart_time, const int first_step)
        {
            m_data_path = data_path;
            m_dt = dt;
            m_cur_time = cur_time;
            m_restart_time = restart_time;
            m_first_step = first_step;
        }

        amrex::Vector<amrex::Real> print_moment(amrex::MultiFab &field, const amrex::Vector<std::string> names,  
                                 const amrex::Vector<int> &moment_orders, SmallDataIO &file, 
                                 const int is_first_step);

    private:
        InflationConfig m_params;
        const int print_mode_functions = 0;
        std::string m_data_path;
        amrex::Real m_dt;
        amrex::Real m_cur_time;
        int m_restart_time;
        int m_first_step;
        const amrex::Real norm;

        std::string make_subdirectory(const std::string base, const std::string dir, 
                                      const int is_first_step) const;

        void assign_statistics_data(amrex::Vector<std::string> &header_storage, const std::string name, 
                                    amrex::Vector<amrex::Real> &data_storage, const amrex::Vector<amrex::Real> data, 
                                    const int component, const int num_comps,
                                    const amrex::Vector<int>::const_iterator itr, 
                                    const amrex::Vector<int>::const_iterator start, 
                                    const int is_first_step);

        void print_power_spectrum(const amrex::cMultiFab &field_array, SmallDataIO &power_spec_file, 
                                  const int component);
        amrex::Real calculate_field_moment_x(const amrex::MultiFab &field, const amrex::Vector<amrex::Real> mean, 
                                      const int moment, const int component);

        void extract_hs_and_R(amrex::MultiFab &hs, amrex::MultiFab &R, 
                              const amrex::MultiFab &state, const bool print_spec);
};

#include "InflationExtraction.impl.hpp"

#endif /* INFLATIONEXTRACTION_HPP_ */
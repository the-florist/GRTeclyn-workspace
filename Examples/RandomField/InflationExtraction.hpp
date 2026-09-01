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
#include <AMReX_GpuContainers.H>

class InflationExtraction
{
    public:
        // Names of diagnostic variables
        static inline const amrex::Vector<std::string> var_names{"R", "hplus", "hcross"};

        // Constructor used in extraction of diagnostics
        InflationExtraction()
        {
            inflt_methods.fill_params();
        }

        static void set_up(int a_state_index);

        // Derive callback (amrex::DeriveFuncMF) that fills plotfile output
        static void compute_mf(amrex::MultiFab &out_mf, int dcomp, int ncomp,
                               const amrex::MultiFab &src_mf,
                               const amrex::Geometry &geomdata,
                               amrex::Real time, const int *bcrec, int level);

        // Main routines
        void derive(const amrex::MultiFab &source,
                    amrex::MultiFab &out, const int dcomp);
        void extract(const amrex::MultiFab &state);

        // Set-up function when printing to data files
        void set_print_params(const std::string data_path,
                              const amrex::Real a_time, const amrex::Real a_dt,
                              const amrex::Real a_restart_time)
        {
            m_data_path  = data_path;
            time         = a_time;
            dt           = a_dt;
            restart_time = a_restart_time;
            first_step   = (time <= dt);
        }

        amrex::Vector<amrex::Real> 
        print_moment(amrex::MultiFab &field, 
                     const amrex::Vector<std::string> names,  
                     const amrex::Vector<int> &moment_orders, 
                     SmallDataIO &file, 
                     const int is_first_step);

    private:
        InflationConfig inflt_methods;
        std::string m_data_path;
        amrex::Real time;
        amrex::Real dt;
        amrex::Real restart_time;
        bool first_step;

        std::string make_subdirectory(const std::string base, const std::string dir, 
                                      const int is_first_step) const;

        void assign_statistics_data(amrex::Vector<std::string> &header_storage,
                                    const std::string name,
                                    amrex::Vector<amrex::Real> &data_storage,
                                    const amrex::Vector<amrex::Real> data,
                                    const int component, const int num_comps,
                                    const amrex::Vector<int>::const_iterator itr, 
                                    const amrex::Vector<int>::const_iterator start, 
                                    const int is_first_step);

        void print_power_spectrum(const amrex::cMultiFab &field_array,
                                  SmallDataIO &power_spec_file, const int component);
                                  
        amrex::Real calculate_field_moment_x(const amrex::MultiFab &field,
                                      const amrex::Vector<amrex::Real> mean,
                                      const int moment, const int component);

        void extract_hs_and_R(amrex::MultiFab &hs, amrex::MultiFab &R, 
                              const amrex::MultiFab &state, const bool print_spec);
};

#include "InflationExtraction.impl.hpp"

#endif /* INFLATIONEXTRACTION_HPP_ */
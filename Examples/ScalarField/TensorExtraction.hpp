/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef TENSOREXTRACTION_HPP_
#define TENSOREXTRACTION_HPP_

#include "RandomField.hpp"

using namespace amrex;

class TensorExtraction : public RandomField
{
    public:
        struct params_t : public RandomField::params_t
        {
            int calc_binned_power_spectrum = 0;   //!< Choose whether to extract the binned power spectrum
            int bin_number = 0;          //!< How many bins to use (capped at N/2)
            int calc_higher_order_statistics = 0; //!< Choose whether to print higher-order statistics on the fields
            int num_orders;
            Vector<int> orders;                   //!< Moment orders to print for extracted fields
        };

        TensorExtraction(params_t a_params, 
                          RandomField::params_t a_random_field_params,
                          InitialBackgroundData::params_t a_background_params) 
                          : m_params(a_params), 
                          RandomField(a_random_field_params, a_background_params)
        {
        }

        void derive(const MultiFab &source, MultiFab &out, int dcomp);
        void extract(const MultiFab &state, const std::string data_path, const Real dt,  
                     const Real cur_time, const int restart_time, const int first_step, const int plot_int);

        void print_tensor_moment(MultiFab &field, const Vector<std::string> names,  
                                 const Vector<int> &moment_orders, SmallDataIO &file, 
                                 const int is_first_step);

    private:
        void print_power_spectrum(cMultiFab &field_array, SmallDataIO &power_spec_file, const int component);
        Real find_field_moment_x(MultiFab &field, const Vector<Real> mean, 
                                 const int moment, const int component);

    protected:
        params_t m_params;
};

#include "TensorExtraction.impl.hpp"

#endif
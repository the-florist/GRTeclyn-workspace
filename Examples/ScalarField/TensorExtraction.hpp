/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef TENSOREXTRACTION_HPP_
#define TENSOREXTRACTION_HPP_

#include "RandomField.hpp"
#include "MatterWeyl4.hpp"

using namespace amrex;

template <class matter_t>
class TensorExtraction : public RandomField, public MatterWeyl4<matter_t>
{
    public:
        using Vars = CCZ4Vars::VarsWithGauge<Real>;
        using Diff2Vars = ADMConformalVars::Diff2VarsNoGauge<Real>;

        TensorExtraction( RandomField::params_t a_params,
                          InitialBackgroundData::params_t a_background_params,
                          matter_t a_matter,
                          const std::array<double, AMREX_SPACEDIM> a_center,
                          const double a_dx,
                          const double a_G_Newton) 
                          : m_params(a_params), m_deriv(a_dx), RandomField(a_params, a_background_params),
                          MatterWeyl4<matter_t>(a_matter, a_center, a_dx, 0, CCZ4RHS<>::USE_CCZ4, a_G_Newton)
        {
        }

        void derive(const MultiFab &source, MultiFab &out, int dcomp);
        void extract_from_Weyl4(const MultiFab &state, const std::string data_path, const Real dt,  
                                 const Real cur_time, const int restart_time, const int first_step, const int plot_int);
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
        RandomField::params_t m_params;
        FourthOrderDerivatives m_deriv; //!< for calculating derivs of vars

};

#include "TensorExtraction.impl.hpp"

#endif
/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef RANDOMFIELD_HPP_
#define RANDOMFIELD_HPP_

#include "Cell.hpp"
#include "InitialScalarData.hpp"
#include "VarsTools.hpp"
#include "FilesystemTools.hpp"
#include <fstream>
#include <random>

#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>
#include <AMReX_FFT.H>
#include <AMReX_Random.H>
#include <AMReX_Print.H>
#include <AMReX_Vector.H>
#include <AMReX_Array.H>

using namespace amrex;

//! Class to create a Gaussian random field, 
//! originally created for 2 massless tensor polarisation fields
//! but will be extended to N IID fields with given masses.
class RandomField 
{
    public:
        // Names of diagnostic variables
        static inline const Vector<std::string> var_names = {"hplus", "hcross"};

        //! A structure for storing parameters essential to this class
        struct params_t 
        {
            int num_scalar_fields;      //!< Number of fields to generate
            int calc_tensor_field;      //!< Determines whether tensor perturbations are calculated
            int use_rand = 1;           //!< Choose whether to use random initial conditions
            int random_seed = 3539263;  //!< Seed for random number generator

            double A;                   //!< Amplitude factor (for basic tests)
            double Mp = 1.;             //!< Energy scale of the problem

            double L_readin;
            int N_readin;               //!< used to read in the private N variable
            int N_fine;                 //!< Fine resolution to downsample from, 
                                        //!< used for convergence testing

            int use_window = 0;         //!< Choose whether to use window function
            double kstar;               //!< window's cut-off mode, measured in units of 2pi/L
            double Delta;               //!< window's width, measured like L/Delta
        };

        RandomField(params_t a_params, InitialBackgroundData::params_t a_background_params)
                : m_params(a_params), m_background_params(a_background_params)
        {
            // Set protected class parameters
            N = m_params.N_readin;
            L = m_params.L_readin;
            norm = m_params.A * pow(2. * M_PI/L, 3.); // Physical FFT normalisation
            tolerance = 1.e-15; // Numerical tolerance, for tests

            // Look-up table 
            // Used to construct polarisation basis tensors
            lut[0][0] = 0;
            lut[0][1] = 1;
            lut[0][2] = 2;
            lut[1][0] = 1;
            lut[1][1] = 3;
            lut[1][2] = 4;
            lut[2][0] = 2;
            lut[2][1] = 4;
            lut[2][2] = 5;
        }

        void init(amrex::MultiFab &state);
        
    protected:
        int N;
        double L;
        int lut[3][3];
        double norm;
        double tolerance;

        // Small functions
        int flip_index(const int indx);
        int invert_index(const int indx);
        int invert_index_with_sign(const int indx);
        void apply_nyquist_conditions(cMultiFab &field);
        Vector<Real> calculate_basis_vector(const IntVect iv, const int which_vector);

        bool is_ghost_index(const IntVect vector);
        std::string make_subdirectory(const std::string base, const std::string dir, const int is_first_step);
        void assign_statistics_data(Vector<std::string> &header_storage, const std::string name, 
                                    Vector<Real> &data_storage, const Vector<Real> data, const int component, const int num_comps,
                                    const Vector<int>::const_iterator itr, const Vector<int>::const_iterator start, 
                                    const int is_first_step);

        // Internal tests
        void Test_is_trace_free(MultiFab &field);

    private:
        // Initialisation sub-routines 
        void make_random_draws(auto &rand_fab, Box &domain);
        GpuComplex<Real> calculate_mode_function(const double km, const std::string spec_type);
        GpuComplex<Real> calculate_random_field(const IntVect iv, const std::string spectrum_type, 
                                                const Real rand_amp, const Real rand_phase);
        GpuComplex<Real> calculate_tensor_initial_conditions(const IntVect iv, const int l, const int p, 
                                                             const GpuComplex<Real> plus_field, 
                                                             const GpuComplex<Real> cross_field);

    //protected:
        const params_t m_params;
        const InitialBackgroundData::params_t m_background_params;
        const std::string m_spec_type;
};

#include "RandomField.impl.hpp"

#endif /* RANDOMFIELD_HPP_ */

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
#include "Potential.hpp"
#include "Tensor.hpp"
#include <fstream>
#include <random>

#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>
#include <AMReX_FFT.H>
#include <AMReX_Random.H>
#include <AMReX_Print.H>
#include <AMReX_Vector.H>
#include <AMReX_Array.H>
#include <AMReX_Reduce.H>
#include <AMReX_ParallelDescriptor.H>

using namespace amrex;

//! Class to create a Gaussian random field, 
//! originally created for 2 massless tensor polarisation fields
//! but will be extended to N IID fields with given masses.
class RandomField 
{
    public:
        // Names of diagnostic variables
        static inline const Vector<std::string> var_names = {"hplus", "hcross", "R"};

        //! A structure for storing parameters essential to this class
        struct params_t 
        {
            // Basic initialisation flags
            int read_from_stoiic = 0;   //!< Whether to read spectrum from stoiic dparams.txt input
            int tensor_init = 0;        //!< Determines whether tensor perturbations are calculated
            int scalar_init = 0;  //!< Read in perturbations from STOIIC dparams
            int use_rand = 1;           //!< Choose whether to use random initial conditions
            int plot_int;
            int window_in_extraction = 0;

            // Grid parameters
            Real L;                   //!< Length of the box
            Real A = 1.;                   //!< Amplitude factor (for basic tests)
            Real Mp = 1.;             //!< Energy scale of the problem
            Real alpha = 0.;          //!< Internal rotation angle in the +/x decomposition basis
            int N_readin;               //!< used to read in the private N variable
            int N_coarse = 0;           //!< Coarse resolution used to set the 
                                        //!< kstar in convergence tests
            int N_fine = 0;            //!< Fine resolution to downsample from, 
                                        //!< used for convergence testing

            // Initial condition options
            int random_seed = 3539263;  //!< Seed for random number generator
            int use_window = 0;         //!< Choose whether to use window function
            Real kstar;               //!< window's cut-off mode, measured in units of 2pi/L
            Real Delta;               //!< window's width, measured like L/Delta

            // Extraction parameters
            int calc_binned_power_spectrum = 0;   //!< Choose whether to extract the binned power spectrum
            int bin_number;                       //!< How many bins to use (capped at N/2)
            int calc_higher_order_statistics = 0; //!< Choose whether to print higher-order statistics on the fields
            int num_orders;                       //!< Number of moments to print (required by vector read-in)
            Vector<int> orders;                   //!< Moment orders to print for extracted fields
            int print_mode_functions = 0;

            // STOIIC read-in structures
            Vector<Real> init_k;                  //!< ks printed by STOIIC, at which Fourier-space fields are provided
            Vector<Vector<Real>> scalar_ps;       //!< Structure: four fields * two components, power spec values
            Vector<Vector<Real>> tensor_ps;       //!< Structure: two fields * two components, power spec values
        };

        // Constructor used when initialising stochastic fields
        RandomField(params_t a_params, InitialBackgroundData::params_t bkgd_params, const Potential::params_t potential_params)
                : m_params(a_params)
        {
            // Compute background potential
            Real V, dV;
            Potential potential(potential_params);
            switch (potential_params.type)
            {
                case 1:
                    potential.quadratic(V, dV, bkgd_params.phi0);
                    break;
                case 4:
                    potential.quadratic_bump(V, dV, bkgd_params.phi0);
                    break;
                case 9:
                    potential.monodromy(V, dV, bkgd_params.phi0);
                    break;
                case 8:
                    potential.USR(V, dV, bkgd_params.phi0);
                    break;
                case 10:
                    potential.punctuated(V, dV, bkgd_params.phi0);
                    break;
                default:
                    Error("RandomField::RandomField, requested " 
                          "potential type is not implemented.");
            }

            // Compute initial Hubble parameter
            H0 = sqrt((8. * M_PI * bkgd_params.G_Newton/3.)*(0.5*pow(bkgd_params.Pi0, 2.) + V));
            phi0 = bkgd_params.phi0;
            pi0 = bkgd_params.Pi0;

            // Set protected class parameters
            N = m_params.N_readin;
            norm = pow(sqrt(2. * M_PI) / m_params.L, 3.); // Physical FFT normalisation
            tolerance = 1.e-10; // Numerical tolerance, for tests

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

        // Constructor used in extraction of diagnostics
        RandomField(params_t a_params)
                : m_params(a_params)
        {
            // Set protected class parameters
            N = m_params.N_readin;
            norm = pow(sqrt(2. * M_PI) / m_params.L, 3.); // Physical FFT normalisation
            tolerance = 1.e-10; // Numerical tolerance, for tests

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
        void derive(const MultiFab &source, MultiFab &out, int dcomp);
        void extract(const MultiFab &state, const std::string data_path, const Real dt,
                     const Real cur_time, const Real restart_time, const int first_step);

        Vector<Real> print_moment(MultiFab &field, const Vector<std::string> names,
                                 const Vector<int> &moment_orders, SmallDataIO &file,
                                 const int is_first_step);

        // Calculates and prints the binned power spectrum of a single
        // real-space component of a MultiFab (e.g. a constraint field
        // populated by a call to derive())
        void print_power_spectrum_of_constraints(MultiFab &field, const int comp, const std::string field_name,
                                           const std::string data_path, const Real dt, const Real cur_time,
                                           const Real restart_time, const int first_step);

    private:
        int N = 0;
        Real H0 = 0.;
        int lut[3][3];
        Real norm = 0.;
        Real tolerance = 0.;
        Real pi0 = 0.;
        Real phi0 = 0.;
        int fft_norm_ft;
        int fft_norm_ift;

        // Small functions
        int flip_index(const int indx, const int m_N);
        int invert_index(const int indx, const int m_N);
        int invert_index_with_sign(const int indx, const int m_N);
        Real get_kmag(IntVect iv, const int m_N);
        Real window_function(const Real kmag);
        Real find_precision_loss(MultiFab &field, int comp, Real bkgd);

        std::string make_subdirectory(const std::string base, const std::string dir, const int is_first_step);
        void assign_statistics_data(Vector<std::string> &header_storage, const std::string name, 
                                    Vector<Real> &data_storage, const Vector<Real> data, const int component, const int num_comps,
                                    const Vector<int>::const_iterator itr, const Vector<int>::const_iterator start, 
                                    const int is_first_step);

        // Tests
        void Test_is_trace_free(MultiFab &field);
        void Test_vector_orthonorm(const IntVect iv, const Vector<Real> mhat, 
                                                                 const Vector<Real> nhat);
        void Test_polarisation_tensor_orthonorm(const IntVect iv, const Tensor<2, Real> eplus,
                                                const Tensor<2, Real> ecross);
        Real calculate_total_power(const cMultiFab& fk, const int m_N);
        void Test_Parsevals_thm(const MultiFab &hx, const cMultiFab &hk, const int m_N);

        // Initialisation routines 
        GpuComplex<Real> calculate_mode_function(const Real km, const int spec_indx);
        GpuComplex<Real> find_in_stoiic(const Real km, const int field_indx, std::string field_type);
        GpuComplex<Real> calculate_random_field(const IntVect iv, const int field_index, 
                                                const Real rand_amp, const Real rand_phase, 
                                                std::string field_type, const int m_N);
        Vector<Real> calculate_basis_vector(const IntVect iv, const int which_vector, const int m_N);
        void apply_nyquist_conditions(cMultiFab &field, const int m_N);
        
        // Extraction routines
        void print_power_spectrum(cMultiFab &field_array, SmallDataIO &power_spec_file, const int component);
        Real find_field_moment_x(MultiFab &field, const Vector<Real> mean, 
                                 const int moment, const int component);

    protected:
        const params_t m_params;
};

#include "RandomField.impl.hpp"

#endif /* RANDOMFIELD_HPP_ */

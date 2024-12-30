/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef RANDOMFIELD_HPP_
#define RANDOMFIELD_HPP_

#include "Cell.hpp"
#include "InitialScalarData.hpp"
#include "VarsTools.hpp"
#include <fstream>

#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>
#include <AMReX_FFT.H>
#include <AMReX_Random.H>

//! Class to create a Gaussian random field, 
//! originally created for 2 massless tensor polarisation fields
//! but will be extended to N IID fields with given masses.
class RandomField 
{
    public:
        //! A structure for storing parameters essential to this class
        struct params_t 
        {
            int num_scalar_fields; //!< Number of fields to generate
            int calc_tensor_field; //!< Determines whether tensor perts are calculated

            double L;          //!< Length of the box
            double A;          //!< Amplitude factor (for basic tests)
            int which_seed;    //!< Which random seed will be chosen (defunct?)

            int N_readin;      //!< used to read in the private N variable
            int N_fine;            //!< Fine resolution to downsample from, 
                                //! used for convergence testing
        };

        RandomField(params_t a_params, InitialBackgroundData::params_t a_background_params, 
                     std::string a_spec_type)
                : m_params(a_params), m_background_params(a_background_params), 
                  m_spec_type(a_spec_type)
        {
        }

        void init_random_field() const;
        
    private:
        int N;                 //<! Grid resolution

    protected:
        const params_t m_params;
        const InitialBackgroundData::params_t m_background_params;
        const std::string m_spec_type;
};


#endif /* RANDOMFIELD_HPP_ */
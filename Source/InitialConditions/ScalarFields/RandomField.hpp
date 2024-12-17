/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef RANDOMFIELD_HPP_
#define RANDOMFIELD_HPP_

#include "Cell.hpp"
#include "Tensor.hpp"
#include "UserVariables.hpp"
#include "InitialScalarData.hpp"
#include "Coordinates.hpp"
#include "VarsTools.hpp"
//#include <random> // needed for random number generator
#include <fstream>

//! Class to create a Gaussian random field, 
//! originally created for 2 massless tensor polarisation fields
//! but will be extended to N IID fields with given masses.
class RandomField 
{
    public:
        //! A structure for storing parameters essential to this class
        struct params_t 
        {
            double num_fields; //!< Number of fields to generate
            double L;          //!< Length of the box
            double A;          //!< Amplitude factor (for basic tests)
            int which_seed;    //!< Which random seed will be chosen (defunct?)
            int Nf;            //!< Fine resolution to downsample from, 
                                //! used for convergence testing
        };

        RandomField(params_t a_params, InitialScalarData::params_t a_background_params, 
                     std::string a_spec_type)
                : m_params(a_params), m_background_params(a_background_params), 
                  m_spec_type(a_spec_type)
        {
        }

        template <class data_t>
        void init_random_field() const;
        
    private:
        int N;                 //<! Grid resolution

    protected:
        const params_t m_params;
        const InitialScalarData::params_t m_background_params;
        const std::string m_spec_type;
}


#endif /* RANDOMFIELD_HPP_ */
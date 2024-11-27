/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef INITIALBACKGROUNDDATA_HPP_
#define INITIALBACKGROUNDDATA_HPP_

#include "Cell.hpp"
#include "Coordinates.hpp"
#include "MatterCCZ4RHS.hpp"
#include "ScalarField.hpp"
#include "StateVariables.hpp" //This files needs NUM_VARS - total no. components
#include "Tensor.hpp"
#include "VarsTools.hpp"
#include "simd.hpp"

//#include "MayDay.H"
//#include <fstream>

class InitialBackgroundData
{
	public:
		struct params_t
		{
			double phi0; //!< Amplitude of k=0 mode of initial SF
			double Pi0;  //!< Amplitude of initial SF velocity
			double m;    //!< SF mass
			double E;    //!< Energy scale [Mp]
		};

		InitialBackgroundData(params_t a_params)
			: m_params(a_params)
		{
		}

		template <class data_t> 
		AMREX_GPU_DEVICE AMREX_FORCE_INLINE void 
		compute(int i, int j, int k, const amrex::Array4<data_t> &cell) const
		{
			const double Mp = 1./m_params.E;

			cell(i, j, k, c_phi) = m_params.phi0;
			cell(i, j, k, c_Pi) = m_params.Pi0;
		}

	protected:
		const params_t m_params;

};


#endif /* INITIALBACKGROUNDDATA_HPP_ */

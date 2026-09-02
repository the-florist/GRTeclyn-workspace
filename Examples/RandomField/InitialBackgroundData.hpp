/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef INITIALBACKGROUNDDATA_HPP_
#define INITIALBACKGROUNDDATA_HPP_

#include "Coordinates.hpp"
#include "StateVariables.hpp"

#include <AMReX_Math.H>

// template <class potential_t>
class InitialBackgroundData
{
	public:
		struct params_t
		{
			amrex::Real phi0; //!< Amplitude of k=0 mode of initial SF
			amrex::Real Pi0;  //!< Amplitude of initial SF velocity
			amrex::Real G_Newton; 

			void fill_params()
			{
				GRParmParse pp;
				pp.get("G_Newton", G_Newton);
				pp.get("init.background_phi", phi0);
				pp.get("init.background_dphi", Pi0);
			}
		};

		InitialBackgroundData() { m_params.fill_params(); }

		template <class data_t> 
		AMREX_GPU_DEVICE AMREX_FORCE_INLINE void 
		compute(int i, int j, int k, const amrex::Array4<data_t> &cell) const
		{
			const amrex::CellData<data_t> &state_cell = cell.cellData(i, j, k);

			// The caller zero-initialises the state, so only the non-zero
			// components of the flat, unit-lapse background are set here.
			state_cell[c_chi]   = 1.0;
			state_cell[c_lapse] = 1.0;
			state_cell[c_h11]   = 1.0;
			state_cell[c_h22]   = 1.0;
			state_cell[c_h33]   = 1.0;

			state_cell[c_phi] = m_params.phi0;
			state_cell[c_Pi]  = m_params.Pi0;

			data_t V, dV;
			m_potential.compute_background_potential(V, dV, m_params.phi0);

			amrex::Real H0 =
				sqrt((8. * amrex::Math::pi<amrex::Real>() * m_params.G_Newton / 3.)
				     * (0.5 * pow(m_params.Pi0, 2.) + V));
			state_cell[c_K] = -3. * H0;
		}
	protected:
		params_t m_params;
		const Potential m_potential;

};


#endif /* INITIALBACKGROUNDDATA_HPP_ */

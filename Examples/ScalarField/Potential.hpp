/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef POTENTIAL_HPP_
#define POTENTIAL_HPP_

#include "simd.hpp"
#include <typeinfo>

class Potential
{
  public:
    struct params_t
    {
		int type;
		amrex::Real param1, param2, param3, param4, param5;

		// Monodromy parameters
		amrex::Real scalar_mass = 0.;
		amrex::Real location;
		amrex::Real width;	
		amrex::Real amplitude;
		amrex::Real period;

		// USR (Prokopec) parameters
		amrex::Real Lambda;
		amrex::Real v = 0.;

		// Punctuated inflation params
		int n;
		amrex::Real lambda;
		// mass

		// Quadratic with Gaussian feature (quadbump) parameters
		// location
		// amp
		// width
    };

  private:
    params_t m_params;

  public:
    //! The constructor
    Potential(params_t a_params) : m_params(a_params) 
	{
		switch (m_params.type)
		{
			case 1:
				m_params.scalar_mass = m_params.param1;
				break;

			case 4:
				m_params.scalar_mass = m_params.param1;
				m_params.location = m_params.param2;
				m_params.amplitude = m_params.param3;
				m_params.width = m_params.param4;
				break;

			case 8:
				m_params.Lambda = m_params.param1;
				m_params.v = m_params.param2;
				break;
			
			case 9:
				m_params.scalar_mass = m_params.param1;
				m_params.location = m_params.param2;
				m_params.width = m_params.param3;
				m_params.amplitude = m_params.param4;
				m_params.period = m_params.param5;
				break;

			case 10:
				m_params.scalar_mass = m_params.param1;
				m_params.n = m_params.param2;
				m_params.lambda = m_params.param3;
				break;

			default:
				amrex::Print() << m_params.type << ", ";
				amrex::Print() << typeid(m_params.type).name() << "\n";
				amrex::Error("Potential::Potential, provided potential "
							 "type is not implemented");
		}
	}

	// Classic quadratic potenital
	template <class data_t>
	AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
	quadratic(data_t &V, data_t &dV, const data_t &phi) const
	{
		if (m_params.scalar_mass == 0)
		{
			amrex::Error("Potential::quadratic, Scalar mass is un-initialised.");
		}

		V = std::pow(m_params.scalar_mass * phi, 2.) / 2.;
		dV = std::pow(m_params.scalar_mass, 2.) * phi;
	}

	// Classic quadratic potenital
	template <class data_t>
	AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
	quadratic_bump(data_t &V, data_t &dV, const data_t &phi) const
	{
		if (m_params.scalar_mass == 0)
		{
			amrex::Error("Potential::quadratic, Scalar mass is un-initialised.");
		}

		amrex::Real feature = m_params.amplitude * std::exp(
							-std::pow((phi - m_params.location) / m_params.width, 2.) 
							/ 2.0);
		
		V = std::pow(m_params.scalar_mass * phi, 2.) * (1.0 + feature) / 2.;
		dV = std::pow(m_params.scalar_mass, 2.) * (phi * (1.0 + feature)
			 - (std::pow(phi, 2.0) * (phi - m_params.location) * feature 
			   / std::pow(m_params.width, 2.0) / 2.0));
	}

	// Monodromy potential, as used in STOIIC and also in arXiv:2403.12811
	template <class data_t>
	AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
	monodromy(data_t &V, data_t &dV, const data_t &phi) const
	{
		if (m_params.scalar_mass == 0)
		{
			amrex::Error("Potential::monodromy, Scalar mass is un-initialised.");
		}

		// Calculate V
		amrex::Real argument = (phi - m_params.location)/m_params.period;
		amrex::Real displaced_argument = (m_params.location - phi + m_params.width)/m_params.period;
		
		amrex::Real envelope = 0.25 * (1. + tanh(argument)) * (1. + tanh(displaced_argument));
		amrex::Real oscillation = cos(argument) - 1.; 	

		V = 0.5 * pow(m_params.scalar_mass * phi, 2.0);
		V += m_params.amplitude * (envelope * oscillation);

		// Calculate dV
		amrex::Real d_envelope = 0.25/m_params.period * 
					((1. + tanh(argument)) * (std::pow(tanh(displaced_argument), 2.) - 1.)
					+ (1. + tanh(displaced_argument)) * (1. - std::pow(tanh(argument), 2.)));
		amrex::Real d_oscillation = -sin(argument)/m_params.period;

		dV = pow(m_params.scalar_mass, 2.0) * phi;
		dV += m_params.amplitude * (envelope * d_oscillation + d_envelope * oscillation);
	}

	// Prokopec USR model, from arXiv:2507.04114
	template <class data_t>
	AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
	USR(data_t &V, data_t &dV, const data_t &phi) const
	{
		if (m_params.v == 0)
		{
			amrex::Error("Potential::USR, USR parameter v is un-initialised.");
		}

		// Calculate V
		amrex::Real fraction = (3. * pow(phi, 2.) 
						 + 2. * sqrt(2.) * phi * m_params.v 
						 + 6. * pow(m_params.v, 2.));

		fraction /= pow(3. * pow(phi, 2.) + 2. * pow(m_params.v, 2.), 2.);
		V = m_params.Lambda * pow(m_params.v, 4.) * pow(phi, 2.) * fraction / 3.;

		// Calculate dV
		fraction = ((2. * m_params.v + sqrt(2.) * phi) 
				  * (pow(phi, 2.) - 2. * pow(m_params.v, 2.)));
		fraction /= pow(2. * pow(m_params.v, 2.) + 3. * pow(phi, 2.), 3.);
		dV = - 2. * m_params.Lambda * pow(m_params.v, 5.) * phi * fraction;
	}

	// Prokopec USR model, from arxiv:0809.3915v2
	template <class data_t>
	AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
	punctuated(data_t &V, data_t &dV, const data_t &phi) const
	{
		if (m_params.scalar_mass == 0 && m_params.lambda == 0)
		{
			amrex::Error("Potential::punctuated, punctuated inflation parameters uninitialised.");
		}

		// Calculate V
		V = pow(m_params.scalar_mass * phi, 2.) / 2.;
		V += m_params.lambda * pow(phi, 2 * (m_params.n - 1)) / 4.;
		V -= sqrt(2. * m_params.lambda * (m_params.n - 1)) * m_params.scalar_mass * pow(phi, m_params.n) / m_params.n;

		// Calculate dV
		dV = pow(m_params.scalar_mass, 2.) * phi;
		dV += m_params.lambda * (m_params.n - 1) * pow(phi, 2 * m_params.n - 3) / 2.;
		dV -= sqrt(2. * m_params.lambda * (m_params.n - 1)) * m_params.scalar_mass * pow(phi, m_params.n - 1);
	}

    //! Set the potential function for the scalar field here
    template <class data_t, template <typename> class vars_t>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
    compute_potential(data_t &V_of_phi, data_t &dVdphi,
                      const vars_t<data_t> &vars) const
    {
		switch (m_params.type)
		{
			case 1:
				quadratic<data_t>(V_of_phi, dVdphi, vars.phi);
				break;
			
			case 4:
				quadratic_bump<data_t>(V_of_phi, dVdphi, vars.phi);
				break;

			case 8:
				USR<data_t>(V_of_phi, dVdphi, vars.phi);
				break;

			case 9:
				monodromy<data_t>(V_of_phi, dVdphi, vars.phi);
				break;

			case 10:
				punctuated<data_t>(V_of_phi, dVdphi, vars.phi);
				break;
			
			default:
				amrex::Print() << m_params.type << "\n";
				amrex::Error("Potential::compute_potential, "
							"requested potential type is not supported.");
		}
    
		/* amrex::Print().SetPrecision(15) << "V: " << V_of_phi << "\n";
		amrex::Print().SetPrecision(15) << "dV: " << dVdphi << "\n";
		amrex::Error(); */
    }
};

#endif /* POTENTIAL_HPP_ */

/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef POTENTIAL_HPP_
#define POTENTIAL_HPP_

#include "GRParmParse.hpp"
#include "ScalarFieldVars.hpp"
#include <typeinfo>

class Potential
{
  public:
    struct params_t
    {
		int type;

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

		void fill_params()
		{
			GRParmParse potential_pp("potential");

			potential_pp.get("type", type);
			if (type != 8) { potential_pp.get("param_1", scalar_mass); }
			switch (type)
			{
				// quadratic
				case 1:
					break;

				// quad+bump
				case 4:
					potential_pp.get("param_2", location);
					potential_pp.get("param_3", amplitude);
					potential_pp.get("param_4", width);
					break;
				
				// USR (Prokopec)
				case 8:
					potential_pp.get("param_1", Lambda);
					potential_pp.get("param_2", v);
					break;

				// Monodromy
				case 9:
					potential_pp.get("param_2", location);
					potential_pp.get("param_3", width);
					potential_pp.get("param_4", amplitude);
					potential_pp.get("param_5", period);
					break;

				// Punctuated inflation
				case 10:
					potential_pp.get("param_2", n);
					potential_pp.get("param_3", lambda);
					break;

				// Inverted quad+bump
				case 11:
					potential_pp.get("param_2", location);
					potential_pp.get("param_3", amplitude);
					potential_pp.get("param_4", width);
					break;

				// quad+step
				case 12:
					potential_pp.get("param_2", amplitude);
					potential_pp.get("param_3", location);
					potential_pp.get("param_4", width);
					break;

				default:
					amrex::Print() << type << ", ";
					amrex::Print() << typeid(type).name() << "\n";
					amrex::Error("Potential::Potential, provided potential "
								"type is not implemented");
			}
		}

		void check_params() const
		{
			// Prokopec (no mass parameter)
			if (type == 8)
			{
				AMREX_ALWAYS_ASSERT_WITH_MESSAGE(v != 0,
					"Potential::USR, USR parameter v is un-initialised");
			}
			// Punctuated needs an extra check
			else if (type == 10)
			{
				AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
					scalar_mass != 0 || lambda != 0,
					"Potential::punctuated, punctuated inflation parameters uninitialised");
			}
			else 
			{
				AMREX_ALWAYS_ASSERT_WITH_MESSAGE(scalar_mass != 0,
					"Potential::quadratic, scalar mass is un-initialised");
			}
						
		}
    };

  private:
	params_t m_params{};

  public:
    //! The constructor
    Potential() { m_params.fill_params(); m_params.check_params(); }

	// Classic quadratic potenital
	template <class data_t>
	AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
	quadratic(data_t &V, data_t &dV, const data_t &phi) const
	{
		V = std::pow(m_params.scalar_mass * phi, 2.) / 2.;
		dV = std::pow(m_params.scalar_mass, 2.) * phi;
	}

	// Classic quadratic potenital
	template <class data_t>
	AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
	quadratic_bump(data_t &V, data_t &dV, const data_t &phi) const
	{
		amrex::Real feature = m_params.amplitude * std::exp(
							-std::pow((phi - m_params.location) / m_params.width, 2.) 
							/ 2.0);
		
		V = std::pow(m_params.scalar_mass * phi, 2.) * (1.0 + feature) / 2.;
		dV = std::pow(m_params.scalar_mass, 2.) * (phi * (1.0 + feature)
			 - (std::pow(phi, 2.0) * (phi - m_params.location) * feature 
			   / std::pow(m_params.width, 2.0) / 4.0));
	}

	// Monodromy potential, as used in STOIIC and also in arXiv:2403.12811
	template <class data_t>
	AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
	monodromy(data_t &V, data_t &dV, const data_t &phi) const
	{
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
	AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
	USR(data_t &V, data_t &dV, const data_t &phi) const
	{
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
	AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
	punctuated(data_t &V, data_t &dV, const data_t &phi) const
	{
		// Calculate V
		V = pow(m_params.scalar_mass * phi, 2.) / 2.;
		V += m_params.lambda * pow(phi, 2 * (m_params.n - 1)) / 4.;
		V -= sqrt(2. * m_params.lambda * (m_params.n - 1)) * m_params.scalar_mass * pow(phi, m_params.n) / m_params.n;

		// Calculate dV
		dV = pow(m_params.scalar_mass, 2.) * phi;
		dV += m_params.lambda * (m_params.n - 1) * pow(phi, 2 * m_params.n - 3) / 2.;
		dV -= sqrt(2. * m_params.lambda * (m_params.n - 1)) * m_params.scalar_mass * pow(phi, m_params.n - 1);
	}

	// Classic quadratic potenital
	template <class data_t>
	AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
	inverted_quadratic_bump(data_t &V, data_t &dV, const data_t &phi) const
	{
		amrex::Real feature = m_params.amplitude * std::exp(
							-std::pow((phi - m_params.location) / m_params.width, 2.) 
							/ 2.0);
		
		V = std::pow(m_params.scalar_mass * phi, 2.) * (1.0 - feature) / 2.;
		dV = std::pow(m_params.scalar_mass, 2.) * (phi * (1.0 - feature)
			 + (std::pow(phi, 2.0) * (phi - m_params.location) * feature 
			   / std::pow(m_params.width, 2.0) / 2.0));
	}

	// Quadratic potential modulated by a tanh step, as used in STOIIC
	// V = 0.5 m^2 phi^2 [1 + c tanh((phi - phi_s) / d)]
	template <class data_t>
	AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
	quadratic_step(data_t &V, data_t &dV, const data_t &phi) const
	{
		amrex::Real step = tanh((phi - m_params.location) / m_params.width);
		amrex::Real d_step = (1.0 - std::pow(step, 2.)) / m_params.width;

		V = std::pow(m_params.scalar_mass * phi, 2.)
			* (1.0 + m_params.amplitude * step) / 2.;
		dV = std::pow(m_params.scalar_mass, 2.) * (phi * (1.0 + m_params.amplitude * step)
			 + std::pow(phi, 2.) * m_params.amplitude * d_step / 2.);
	}

    //! Set the potential function for the scalar field here
    template <class data_t, template <typename> class vars_t>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
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

			case 11:
				inverted_quadratic_bump<data_t>(V_of_phi, dVdphi, vars.phi);
				break;

			case 12:
				quadratic_step<data_t>(V_of_phi, dVdphi, vars.phi);
				break;
			
			default:
				AMREX_ALWAYS_ASSERT_WITH_MESSAGE(false,
					"Potential::compute_potential, "
					"requested potential type is not supported");
		}
    }

    //! Concrete overload used by the ScalarField matter class
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
    compute(amrex::Real &V_of_phi, amrex::Real &dVdphi,
                      const ScalarFieldVars &vars) const
    {
		compute_background_potential(V_of_phi, dVdphi, vars.phi());
    }

	//! Set the potential function for the scalar field mean value
    template <class data_t>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
    compute_background_potential(data_t &V_of_phi, data_t &dVdphi,
                      			const data_t phi) const
    {
		switch (m_params.type)
		{
			case 1:
				quadratic<data_t>(V_of_phi, dVdphi, phi);
				break;
			case 4:
				quadratic_bump<data_t>(V_of_phi, dVdphi, phi);
				break;
			case 8:
				USR<data_t>(V_of_phi, dVdphi, phi);
				break;
			case 9:
				monodromy<data_t>(V_of_phi, dVdphi, phi);
				break;
			case 10:
				punctuated<data_t>(V_of_phi, dVdphi, phi);
				break;
			case 11:
				inverted_quadratic_bump<data_t>(V_of_phi, dVdphi, phi);
				break;
			case 12:
				quadratic_step<data_t>(V_of_phi, dVdphi, phi);
				break;
			default:
				AMREX_ALWAYS_ASSERT_WITH_MESSAGE(false,
    				"Potential::compute_background_potential, "
					" potential type not supported");
		}
	}
};

#endif /* POTENTIAL_HPP_ */

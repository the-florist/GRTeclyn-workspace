/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef POTENTIAL_HPP_
#define POTENTIAL_HPP_

#include "GRParmParse.hpp"
#include "ScalarFieldVars.hpp"

#include <AMReX_GpuQualifiers.H>
#include <AMReX_REAL.H>

#include <string>

class Potential
{
  public:
    enum class Type
    {
        Quadratic,
        QuadraticBump,
        USR,
        Monodromy,
        Punctuated,
        InvertedQuadraticBump,
        QuadraticStep
    };

    struct params_t
    {
        Type type{Type::Quadratic}; //!< Which model to evaluate

        amrex::Real scalar_mass{1.0}; //!< Inflaton mass, read from
                                      //!< scalar_field.scalar_mass

        amrex::Real feature_amplitude{0.};
        amrex::Real feature_location{0.};
        amrex::Real feature_width{1.};
        amrex::Real feature_period{1.}; //!< Monodromy only

        amrex::Real usr_lambda{0.};
        amrex::Real usr_v{0.};

        int punctuated_n{2};
        amrex::Real punctuated_lambda{0.};

        static bool type_from_string(const std::string &name, Type &type)
        {
            if (name == "quadratic")
            {
                type = Type::Quadratic;
            }
            else if (name == "quadratic_bump")
            {
                type = Type::QuadraticBump;
            }
            else if (name == "usr")
            {
                type = Type::USR;
            }
            else if (name == "monodromy")
            {
                type = Type::Monodromy;
            }
            else if (name == "punctuated")
            {
                type = Type::Punctuated;
            }
            else if (name == "inverted_quadratic_bump")
            {
                type = Type::InvertedQuadraticBump;
            }
            else if (name == "quadratic_step")
            {
                type = Type::QuadraticStep;
            }
            else
            {
                return false;
            }

            return true;
        }

        static void check_params()
        {
            GRParmParse potential_pp("potential");
            std::string type_name{"quadratic"};
            potential_pp.queryAdd("type", type_name);

            Type type{Type::Quadratic};
            if (!type_from_string(type_name, type))
            {
                potential_pp.error(
                    "type", "must be one of quadratic, quadratic_bump, usr, "
                            "monodromy, punctuated, inverted_quadratic_bump, "
                            "quadratic_step");
            }

            GRParmParse scalar_field_pp("scalar_field");
            amrex::Real scalar_mass{1.0};
            scalar_field_pp.queryAdd("scalar_mass", scalar_mass);
            if (scalar_mass < 0.0)
            {
                scalar_field_pp.error("scalar_mass", "must be >= 0.0");
            }

            // The USR model is the only one whose mass scale is set by
            // something other than scalar_mass.
            if (scalar_mass == 0.0 && type != Type::USR)
            {
                scalar_field_pp.error(
                    "scalar_mass", "must be non-zero for this potential.type");
            }

            check_model_params(type);

            GRParmParse geometry_pp("geometry");
            amrex::Real coarsest_dx{};
            geometry_pp.get("coarsest_dx", coarsest_dx);

            GRParmParse evolution_pp("evolution");
            amrex::Real dt_multiplier{};
            evolution_pp.get("dt_multiplier", dt_multiplier);
            if (scalar_mass >= 0.2 / coarsest_dx / dt_multiplier)
            {
                scalar_field_pp.warning(
                    "scalar_mass",
                    "oscillations of the scalar field may not be resolved on "
                    "the coarsest level");
            }
        }

        //! Validates the parameters belonging to the selected model only, so
        //! that an unused parameter left over in a parameter file cannot fail
        //! the run.
        static void check_model_params(const Type type)
        {
            GRParmParse potential_pp("potential");

            const bool uses_feature_width =
                (type == Type::QuadraticBump ||
                 type == Type::InvertedQuadraticBump ||
                 type == Type::QuadraticStep || type == Type::Monodromy);

            if (uses_feature_width)
            {
                amrex::Real feature_width{1.};
                potential_pp.queryAdd("feature_width", feature_width);
                if (feature_width == 0.)
                {
                    potential_pp.error("feature_width", "must be non-zero");
                }
            }

            if (type == Type::Monodromy)
            {
                amrex::Real feature_period{1.};
                potential_pp.queryAdd("feature_period", feature_period);
                if (feature_period == 0.)
                {
                    potential_pp.error("feature_period", "must be non-zero");
                }
            }

            if (type == Type::USR)
            {
                amrex::Real usr_v{0.};
                potential_pp.queryAdd("usr_v", usr_v);
                if (usr_v == 0.)
                {
                    potential_pp.error("usr_v", "must be non-zero");
                }
            }

            if (type == Type::Punctuated)
            {
                amrex::Real punctuated_lambda{0.};
                potential_pp.queryAdd("punctuated_lambda", punctuated_lambda);
                if (punctuated_lambda == 0.)
                {
                    potential_pp.error("punctuated_lambda", "must be non-zero");
                }

                int punctuated_n{2};
                potential_pp.queryAdd("punctuated_n", punctuated_n);
                if (punctuated_n < 2)
                {
                    potential_pp.error("punctuated_n", "must be >= 2");
                }
            }
        }

        void fill_params()
        {
            GRParmParse scalar_field_pp("scalar_field");
            scalar_field_pp.get("scalar_mass", scalar_mass);

            GRParmParse potential_pp("potential");
            std::string type_name{"quadratic"};
            potential_pp.query("type", type_name);
            if (!type_from_string(type_name, type))
            {
                potential_pp.error("type", "is not a recognised potential");
            }

            potential_pp.query("feature_amplitude", feature_amplitude);
            potential_pp.query("feature_location", feature_location);
            potential_pp.query("feature_width", feature_width);
            potential_pp.query("feature_period", feature_period);

            potential_pp.query("usr_lambda", usr_lambda);
            potential_pp.query("usr_v", usr_v);

            potential_pp.query("punctuated_n", punctuated_n);
            potential_pp.query("punctuated_lambda", punctuated_lambda);
        }
    };

    Potential() { m_params.fill_params(); }

    AMREX_GPU_HOST_DEVICE
    AMREX_FORCE_INLINE explicit Potential(params_t a_params)
        : m_params(a_params)
    {
    }

    //! Evaluates the potential for the evolved scalar field.
    // NOLINTBEGIN(bugprone-easily-swappable-parameters)
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
    compute_potential(amrex::Real &V_of_phi, amrex::Real &dVdphi,
                      const ScalarFieldVars &vars) const
    {
        eval(V_of_phi, dVdphi, vars.phi());
    }
    // NOLINTEND(bugprone-easily-swappable-parameters)

    //! Evaluates the same potential for a homogeneous background value. Used
    //! on the host when setting the initial Friedmann background.
    // NOLINTBEGIN(bugprone-easily-swappable-parameters)
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
    compute_background_potential(amrex::Real &V_of_phi, amrex::Real &dVdphi,
                                 const amrex::Real &phi) const
    {
        eval(V_of_phi, dVdphi, phi);
    }
    // NOLINTEND(bugprone-easily-swappable-parameters)

  private:
    params_t m_params{};

    //! Single evaluation point for every model, so the host background and
    //! the evolved field can never drift apart.
    // NOLINTBEGIN(bugprone-easily-swappable-parameters)
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
    eval(amrex::Real &V_of_phi, amrex::Real &dVdphi,
         const amrex::Real phi) const
    {
        switch (m_params.type)
        {
        case Type::Quadratic:
            quadratic(V_of_phi, dVdphi, phi);
            break;
        case Type::QuadraticBump:
            quadratic_bump(V_of_phi, dVdphi, phi, 1.0);
            break;
        case Type::InvertedQuadraticBump:
            quadratic_bump(V_of_phi, dVdphi, phi, -1.0);
            break;
        case Type::USR:
            usr(V_of_phi, dVdphi, phi);
            break;
        case Type::Monodromy:
            monodromy(V_of_phi, dVdphi, phi);
            break;
        case Type::Punctuated:
            punctuated(V_of_phi, dVdphi, phi);
            break;
        case Type::QuadraticStep:
            quadratic_step(V_of_phi, dVdphi, phi);
            break;
        }
    }
    // NOLINTEND(bugprone-easily-swappable-parameters)

    //! V = (1/2) m^2 phi^2
    // NOLINTBEGIN(bugprone-easily-swappable-parameters)
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
    quadratic(amrex::Real &V_of_phi, amrex::Real &dVdphi,
              const amrex::Real phi) const
    {
        const amrex::Real mass_times_phi = m_params.scalar_mass * phi;
        V_of_phi = 0.5 * mass_times_phi * mass_times_phi;
        dVdphi   = m_params.scalar_mass * m_params.scalar_mass * phi;
    }
    // NOLINTEND(bugprone-easily-swappable-parameters)

    //! Quadratic potential modulated by a Gaussian feature,
    //! V = (1/2) m^2 phi^2 [1 + s A exp(-(phi - phi_0)^2 / 2 w^2)],
    //! where the sign s is +1 for a bump and -1 for a dip.
    // NOLINTBEGIN(bugprone-easily-swappable-parameters)
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
    quadratic_bump(amrex::Real &V_of_phi, amrex::Real &dVdphi,
                   const amrex::Real phi, const amrex::Real sign) const
    {
        const amrex::Real offset =
            (phi - m_params.feature_location) / m_params.feature_width;
        const amrex::Real feature =
            sign * m_params.feature_amplitude * exp(-0.5 * offset * offset);

        const amrex::Real mass_squared =
            m_params.scalar_mass * m_params.scalar_mass;

        // d(feature)/dphi = -feature (phi - phi_0) / w^2, so
        // dV = m^2 [phi (1 + f) - phi^2 f (phi - phi_0) / (2 w^2)]
        const amrex::Real d_feature =
            -feature * (phi - m_params.feature_location) /
            (m_params.feature_width * m_params.feature_width);

        V_of_phi = 0.5 * mass_squared * phi * phi * (1.0 + feature);
        dVdphi   = mass_squared *
                   (phi * (1.0 + feature) + 0.5 * phi * phi * d_feature);
    }
    // NOLINTEND(bugprone-easily-swappable-parameters)

    //! Monodromy potential, as used in STOIIC and in arXiv:2403.12811
    // NOLINTBEGIN(bugprone-easily-swappable-parameters)
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
    monodromy(amrex::Real &V_of_phi, amrex::Real &dVdphi,
              const amrex::Real phi) const
    {
        const amrex::Real period   = m_params.feature_period;
        const amrex::Real argument = (phi - m_params.feature_location) / period;
        const amrex::Real displaced_argument =
            (m_params.feature_location - phi + m_params.feature_width) / period;

        const amrex::Real tanh_arg       = tanh(argument);
        const amrex::Real tanh_displaced = tanh(displaced_argument);

        const amrex::Real envelope =
            0.25 * (1. + tanh_arg) * (1. + tanh_displaced);
        const amrex::Real oscillation = cos(argument) - 1.;

        const amrex::Real d_envelope =
            0.25 / period *
            ((1. + tanh_arg) * (tanh_displaced * tanh_displaced - 1.) +
             (1. + tanh_displaced) * (1. - tanh_arg * tanh_arg));
        const amrex::Real d_oscillation = -sin(argument) / period;

        const amrex::Real mass_squared =
            m_params.scalar_mass * m_params.scalar_mass;

        V_of_phi = 0.5 * mass_squared * phi * phi +
                   m_params.feature_amplitude * envelope * oscillation;
        dVdphi   = mass_squared * phi +
                   m_params.feature_amplitude *
                       (envelope * d_oscillation + d_envelope * oscillation);
    }
    // NOLINTEND(bugprone-easily-swappable-parameters)

    //! Prokopec ultra-slow-roll model, from arXiv:2507.04114
    // NOLINTBEGIN(bugprone-easily-swappable-parameters)
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
    usr(amrex::Real &V_of_phi, amrex::Real &dVdphi, const amrex::Real phi) const
    {
        const amrex::Real v_val       = m_params.usr_v;
        const amrex::Real phi_squared = phi * phi;
        const amrex::Real v_squared   = v_val * v_val;

        const amrex::Real denominator = 3. * phi_squared + 2. * v_squared;

        amrex::Real fraction  = 3. * phi_squared +
                                2. * std::sqrt(2.) * phi * v_val +
                                6. * v_squared;
        fraction             /= denominator * denominator;
        V_of_phi = m_params.usr_lambda * v_squared * v_squared * phi_squared *
                   fraction / 3.;

        fraction =
            (2. * v_val + std::sqrt(2.) * phi) * (phi_squared - 2. * v_squared);
        fraction /= denominator * denominator * denominator;
        dVdphi    = -2. * m_params.usr_lambda * v_squared * v_squared * v_val *
                    phi * fraction;
    }
    // NOLINTEND(bugprone-easily-swappable-parameters)

    //! Punctuated inflation, from arXiv:0809.3915
    // NOLINTBEGIN(bugprone-easily-swappable-parameters)
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
    punctuated(amrex::Real &V_of_phi, amrex::Real &dVdphi,
               const amrex::Real phi) const
    {
        const amrex::Real n_val  = m_params.punctuated_n;
        const amrex::Real lambda = m_params.punctuated_lambda;
        const amrex::Real mass   = m_params.scalar_mass;

        const amrex::Real cross_term =
            std::sqrt(2. * lambda * (n_val - 1.)) * mass;

        V_of_phi = 0.5 * mass * mass * phi * phi +
                   lambda * std::pow(phi, 2. * (n_val - 1.)) / 4. -
                   cross_term * std::pow(phi, n_val) / n_val;

        dVdphi = mass * mass * phi +
                 lambda * (n_val - 1.) * std::pow(phi, 2. * n_val - 3.) / 2. -
                 cross_term * std::pow(phi, n_val - 1.);
    }
    // NOLINTEND(bugprone-easily-swappable-parameters)

    //! Quadratic potential modulated by a tanh step, as used in STOIIC,
    //! V = (1/2) m^2 phi^2 [1 + c tanh((phi - phi_s) / d)]
    // NOLINTBEGIN(bugprone-easily-swappable-parameters)
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
    quadratic_step(amrex::Real &V_of_phi, amrex::Real &dVdphi,
                   const amrex::Real phi) const
    {
        const amrex::Real step =
            tanh((phi - m_params.feature_location) / m_params.feature_width);
        const amrex::Real d_step = (1.0 - step * step) / m_params.feature_width;

        const amrex::Real mass_squared =
            m_params.scalar_mass * m_params.scalar_mass;
        const amrex::Real amplitude = m_params.feature_amplitude;

        V_of_phi = 0.5 * mass_squared * phi * phi * (1.0 + amplitude * step);
        dVdphi   = mass_squared * (phi * (1.0 + amplitude * step) +
                                   0.5 * phi * phi * amplitude * d_step);
    }
    // NOLINTEND(bugprone-easily-swappable-parameters)
};

#endif /* POTENTIAL_HPP_ */

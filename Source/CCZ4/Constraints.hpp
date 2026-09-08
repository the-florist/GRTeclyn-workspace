/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

// This compute class calculates Hamiltonian and Momentum constraints

#ifndef CONSTRAINTS_HPP_
#define CONSTRAINTS_HPP_

// GRTeclyn includes
#include "CCZ4Geometry.hpp"
#include "CCZ4Vars.hpp"
#include "DimensionDefinitions.hpp"
#include "FourthOrderDerivatives.hpp"
#include "GRInterval.hpp"
#include "Interval.hpp"

#include "Tensor.hpp"
#include "TensorAlgebra.hpp"

// AMReX includes
#include <AMReX_MultiFab.H>
#include <AMReX_REAL.H>

// System includes
#include <array>

class Constraints
{
  public:
    static inline const std::string name = "constraints";

    /// Variable names
    static inline const amrex::Vector<std::string> var_names = {"Ham", "Mom1",
                                                                "Mom2", "Mom3"};

    static inline const amrex::Vector<std::string> var_names_norm = {"Ham",
                                                                     "Mom"};

    /// Variable names for the absolute-value calibration terms (Mom stored
    /// as a norm); only valid alongside var_names_norm
    static inline const amrex::Vector<std::string> var_names_abs_terms_norm =
        {"Ham_abs_terms", "Mom_abs_terms"};

    /// Variable names for the relative constraint violation,
    /// i.e. |Ham| / Ham_abs_terms and |Mom| / Mom_abs_terms
    static inline const amrex::Vector<std::string> var_names_relative = {
        "Ham_absrel", "Mom_absrel"};

    /// Struct for Constraints
    struct constraints_t
    {
        amrex::Real Ham{};
        amrex::Real Ham_abs_terms{};
        Tensor::Rank1 Mom{};
        Tensor::Rank1 Mom_abs_terms{};
    };

    // Constructor which allows specifying Ham and Mom vars
    // if the interval of a_c_Moms is of size 1, then
    // sqrt(Mom1^2 + Mom2^2 + Mom3^2) is stored in that variable
    // ...abs_terms stores the absolute value of the individual terms in the
    // conformally decomposed expressions which can be used in to normalize
    // the constraint violations
    // Any zero-length Interval or negative var is not calculated
    // If a positive interval is passed for one of a_c_Moms or
    // a_c_moms_abs_terms then it must have length consistent with
    // s_calc_mom_norm
    AMREX_FORCE_INLINE
    Constraints(amrex::Real dx, int a_c_Ham, const Interval &a_c_Moms,
                int a_c_Ham_abs_terms              = -1,
                const Interval &a_c_Moms_abs_terms = Interval(),
                amrex::Real cosmological_constant  = 0.0);

    AMREX_FORCE_INLINE AMREX_GPU_DEVICE void
    operator()(int ix, int iy, int iz,
               const amrex::Array4<amrex::Real> &constraints,
               const amrex::Array4<amrex::Real const> &state) const;

    // Computes the relative constraint violation |Ham| / Ham_abs_terms and
    // |Mom| / Mom_abs_terms from the values already written by operator()
    // into the constraints array at m_c_Ham/m_c_Ham_abs_terms and the
    // (norm-form) Mom equivalents. Must be called after operator() has run
    // on the same cell. Divisions are floored to avoid NaN/Inf where all
    // constraint terms vanish (e.g. exact vacuum).
    AMREX_FORCE_INLINE AMREX_GPU_DEVICE void
    compute_absrel(int ix, int iy, int iz,
                   const amrex::Array4<amrex::Real> &constraints,
                   int c_Ham_absrel, int c_Mom_absrel) const;

    /// Adds the constraints to the derive list
    /// Call in variableSetUp()
    /// a_calc_abs_terms additionally computes and stores Ham_abs_terms,
    /// Mom_abs_terms, Ham_absrel and Mom_absrel; it requires
    /// a_calc_mom_norm == true, since the relative constraint needs the
    /// Mom-norm layout
    AMREX_FORCE_INLINE static void set_up(int a_state_index,
                                          bool a_calc_mom_norm  = false,
                                          bool a_calc_abs_terms = false);

    // Has signature of DeriveFuncMF so that it can be stored in the derive_lst
    AMREX_FORCE_INLINE static void
    compute_mf(amrex::MultiFab &out_mf, int dcomp, int ncomp,
               const amrex::MultiFab &src_mf, const amrex::Geometry &geomdata,
               amrex::Real /*time*/, const int * /*bcrec*/, int /*level*/);

  protected:
    // NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
    static inline bool s_calc_mom_norm =
        false; // set to true with set_up() to store just sqrt(Mom1^2 + Mom2^2 +
               // Mom3^2) instead of Mom1, Mom2, Mom3 separately
    // NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
    static inline bool s_calc_abs_terms =
        false; // set to true with set_up() to additionally store
               // Ham_abs_terms, Mom_abs_terms, Ham_absrel, Mom_absrel
    FourthOrderDerivatives m_deriv;
    int m_c_Ham;
    Interval m_c_Moms;
    int m_c_Ham_abs_terms = -1;
    Interval m_c_Moms_abs_terms;
    amrex::Real m_cosmological_constant;

    [[nodiscard]]
    AMREX_FORCE_INLINE AMREX_GPU_DEVICE constraints_t constraint_equations(
        const CCZ4Vars &vars, const Tensor::Rank1 &d1_chi,
        const Tensor::Rank2 &d1_Gamma, const Tensor::Sym12Rank3 &d1_h,
        const Tensor::Rank1 &d1_K, const Tensor::Sym12Rank3 &d1_A,
        const Tensor::Sym12Rank2 &d2_chi, const Tensor::Sym12Sym34Rank4 &d2_h,
        const Tensor::Rank2 &h_UU, const chris_t &chris) const;

    AMREX_FORCE_INLINE AMREX_GPU_DEVICE void
    store_vars(const constraints_t &out,
               const amrex::CellData<amrex::Real> &current_cell) const;
};

#include "Constraints.impl.hpp"

#endif /* CONSTRAINTS_HPP_ */

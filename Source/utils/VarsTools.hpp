/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef VARSTOOLS_HPP_
#define VARSTOOLS_HPP_

// Our includes
#include "GRInterval.hpp"
#include "StateVariables.hpp"
#include "Tensor.hpp"

// AMReX includes
#include "AMReX_Print.H" //Gives us amrex::Print()

namespace VarsTools
{
template <typename mapping_function_t, typename data_t>
AMREX_GPU_HOST_DEVICE void
define_enum_mapping(mapping_function_t mapping_function, const int &ivar,
                    data_t &scalar)
{
    mapping_function(ivar, scalar);
}

template <typename mapping_function_t, typename data_t, int start_var,
          int end_var>
AMREX_GPU_HOST_DEVICE void
define_enum_mapping(mapping_function_t mapping_function,
                    const GRInterval<start_var, end_var> interval,
                    Tensor<1, data_t, end_var - start_var + 1> &tensor)
{
    for (int ivar = 0; ivar < interval.size(); ++ivar)
    {
        mapping_function(start_var + ivar, tensor[ivar]);
    }
}

template <typename mapping_function_t, typename data_t, int start_var,
          int end_var>
AMREX_GPU_HOST_DEVICE void
define_symmetric_enum_mapping(mapping_function_t mapping_function,
                              const GRInterval<start_var, end_var> interval,
                              Tensor<2, data_t> &tensor)
{
    static_assert(interval.size() ==
                      DEFAULT_TENSOR_DIM * (DEFAULT_TENSOR_DIM + 1) / 2,
                  "Interval has wrong size");
#if DEFAULT_TENSOR_DIM == 3
    mapping_function(start_var, tensor[0][0]);

    mapping_function(start_var + 1, tensor[0][1]);
    mapping_function(start_var + 1, tensor[1][0]);

    mapping_function(start_var + 2, tensor[0][2]);
    mapping_function(start_var + 2, tensor[2][0]);

    mapping_function(start_var + 3, tensor[1][1]);

    mapping_function(start_var + 4, tensor[1][2]);
    mapping_function(start_var + 4, tensor[2][1]);

    mapping_function(start_var + 5, tensor[2][2]);
#else
#error DEFAULT_TENSOR_DIM not equal to three not implemented yet...
#endif
}

//--> Begin: Helper for the assign function
template <class nested_template> struct strip_nested_template;

template <template <typename> class outermost_layer, class inner_part>
struct strip_nested_template<outermost_layer<inner_part>>
{
    using type = inner_part;
};
//<-- End: Helper for the assign function

/// Writes data directly into all variables
/**if this variables has multiple components (e.g. if it is an array of
 *derivatives) the data can be written directly into these components by
 *specifying an arbitrary number of icomps
 */
template <class vars_t, typename value_t>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void assign(vars_t &vars,
                                                     const value_t &value)
{
    // The template magic below is needed to make sure that we can write
    // assign(vars, 0.)  and 0. gets correctly cast from double to simd<double>
    // if necessary.
    using data_t = typename strip_nested_template<vars_t>::type;
    vars.enum_mapping([&value](const int & /*ivar*/, data_t &var)
                      { var = static_cast<data_t>(value); });
}

/// Prints all elements of the vars element with component names
/// (Very useful for debugging)
template <template <typename> class vars_t, typename data_t>
void print(const vars_t<data_t> &vars)
{
    vars.enum_mapping(
        [](const int &ivar, data_t &var) {
            amrex::Print() << StateVariables::names[ivar] << ": " << var
                           << "\n";
        });
}
} // namespace VarsTools

#endif /* VARSTOOLS_HPP_ */

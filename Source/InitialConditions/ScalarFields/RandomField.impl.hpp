/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */


#if !defined(RANDOMFIELD_HPP_)
#error "This file should only be included via RandomField.hpp"
#endif

#ifndef RANDOMFIELD_IMPL_HPP_
#define RANDOMFIELD_IMPL_HPP_

/****
    GpuParams construction
****/

inline RandomField::GpuParams RandomField::make_gpu_params() const
{
    GpuParams gp;
    gp.N         = N;
    gp.H0        = H0;
    gp.norm      = norm;
    gp.tolerance = tolerance;
    for(int l = 0; l < 3; l++)
        for(int p = 0; p < 3; p++)
            gp.lut[l][p] = lut[l][p];

    gp.scalar_init      = m_params.scalar_init;
    gp.tensor_init      = m_params.tensor_init;
    gp.use_rand         = m_params.use_rand;
    gp.use_window       = m_params.use_window;
    gp.N_coarse         = m_params.N_coarse;
    gp.read_from_stoiic = m_params.read_from_stoiic;
    gp.random_seed      = m_params.random_seed;
    gp.L                = m_params.L;
    gp.A                = m_params.A;
    gp.alpha            = m_params.alpha;
    gp.Delta            = m_params.Delta;
    // Stoiic device pointers are left nullptr; init() fills them when needed.
    return gp;
}

/****
    Index utilities (static, AMREX_GPU_HOST_DEVICE)
****/

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int RandomField::flip_index(const int indx, const int m_N)
{
    return std::abs(m_N - indx);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int RandomField::invert_index(const int indx, const int m_N)
{
    return (int)(m_N/2 - std::abs(m_N/2 - indx));
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int RandomField::invert_index_with_sign(const int indx, const int m_N)
{
    if(indx <= m_N/2) { return indx; }
    else { return std::abs(m_N/2 - indx) - m_N/2; }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real RandomField::get_kmag(IntVect iv, const int m_N, const Real L)
{
    const int i = iv[0];
    const int j = invert_index(iv[1], m_N);
    const int k = invert_index(iv[2], m_N);
    return std::sqrt(static_cast<Real>(i*i + j*j + k*k)) * 2. * M_PI / L;
}

/****
    Small host-only helpers
****/

inline std::string RandomField::make_subdirectory(const std::string base, const std::string dir, const int is_first_step)
{
    std::string new_path = base+"../"+dir+"/";
    if(is_first_step)
    {
        if (FilesystemTools::directory_exists(base)) { FilesystemTools::mkdir_recursive(new_path); }
        else 
        { 
            Print() << "RandomField::make_subdirectory, Directory creation failed for " << new_path << "\n";
            Error("RandomField::extract Data directory has not been created."); 
        }
    }
    return new_path;
}

// Creates a custom data file layout 
inline void RandomField::assign_statistics_data(Vector<std::string> &header_storage, const std::string name, 
                            Vector<Real> &data_storage, const Vector<Real> data, const int component, const int num_comps,
                            const Vector<int>::const_iterator itr, const Vector<int>::const_iterator start, 
                            const int is_first_step)
{
    int loc = component + num_comps*(itr - start);
    if(is_first_step) 
    { 
        header_storage[loc] =  name; 
    }
    data_storage[loc] = data[component];
}

// Written by Gemini, edited by me
inline Real RandomField::find_precision_loss(MultiFab &field, int comp, Real bkgd)
{
    // 1. Initialize the reduction operator for a 'Minimum' operation
    amrex::ReduceOps<amrex::ReduceOpMin> reduce_op;
    amrex::ReduceData<amrex::Real> reduce_data(reduce_op);
    using ReduceTuple = typename decltype(reduce_data)::Type;

    // 2. Loop over the MultiFab (GPU and CPU safe)
    for (amrex::MFIter mfi(field); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.fabbox();
        auto const& arr = field.array(mfi);

        reduce_op.eval(bx, reduce_data,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple {
                // Return the absolute value of the cell
                return { std::abs(arr(i, j, k, comp)) };
            });
    }

    // 3. Extract the local minimum on this MPI rank
    amrex::Real min_abs_val = amrex::get<0>(reduce_data.value());

    // 4. Perform a collective MPI reduction to find the global minimum across all processors
    amrex::ParallelDescriptor::ReduceRealMin(min_abs_val);

    int p_field = std::round(std::log10(min_abs_val));
    int p_bkgd = std::round(std::log10(std::abs(bkgd)));

    if (p_bkgd + p_field > 0)
    {
        Print() << bkgd << ", " << min_abs_val << "\n";
        Print() << p_bkgd << ", " << p_field << "\n";
        Error("RandomField::find_precision_loss, field may be non-perturbative.");
    }

    return pow(10., p_bkgd + p_field);
}

/****
    Tests
****/

// Test that the input tensor field (config space) is trace free (global).
// Uses AMREX_ASSERT inside the kernel; diagnostic prints stripped for GPU compat.
inline void RandomField::Test_is_trace_free(MultiFab &field, const GpuParams& gp)
{
    if (field.nComp() != 6)
    {
        Error("RandomField::Test_is_trace_free, input field is not a tensor field");
    }

    for (MFIter mfi(field); mfi.isValid(); ++mfi) 
    {
        Array4<Real> const& field_ptr = field.array(mfi);
        const Box& bx = mfi.fabbox();

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real sum = 0.;
            for(int l = 0; l < 3; l++)
            {
                sum += field_ptr(i, j, k, gp.lut[l][l]);
            }
            AMREX_ASSERT_WITH_MESSAGE(std::abs(sum) <= gp.tolerance,
                "RandomField::Test_is_trace_free Trace-free test failed.");
        });
    }
}

// Test that the input vectors are orthonormal (host only)
inline void RandomField::Test_vector_orthonorm(const IntVect iv,
                                               const GpuArray<Real, 3>& mhat,
                                               const GpuArray<Real, 3>& nhat)
{
    // Confirm basis vectors are orthonormal
    if (iv != IntVect{0, 0, 0})
    {
        Real dot1 = 0., dot2 = 0., cross1 = 0.;
        for(int l = 0; l < 3; l++)
        {
            dot1   += mhat[l] * mhat[l];
            dot2   += nhat[l] * nhat[l];
            cross1 += mhat[l] * nhat[l];
        }

        if(std::abs(dot1 - 1.) > tolerance || std::abs(dot2 - 1.) > tolerance || cross1 > tolerance)
        {
            Print() << "Location: " << iv << "\n";
            Print() << "Dot products: " << dot1 << ", " << dot2 << "\n";
            Print() << "Cross product: " << cross1 << "\n";
            Print() << "Vectors: \n";
            for(int l = 0; l < 3; l++)
            {
                Print() << l << ", " << mhat[l] << ", " << nhat[l] << "\n";
            }
            Error("RandomField::Test_vector_orthonorm: Basis vectors are not orthonormal here");
        }
    }
}

// Test that the input basis tensors, and their rotated counterparts, are orthonormal (host only)
inline void RandomField::Test_polarisation_tensor_orthonorm(const IntVect iv, const Tensor<2, Real> eplus,
                                                                              const Tensor<2, Real> ecross)
{
    Vector<Real> conds(3, 0.);

    for (int l=0; l<3; l++) for (int p=0; p<3; p++)
    {
        conds[0] += eplus[l][p] * eplus[l][p];
        conds[1] += eplus[l][p] * ecross[l][p];
        conds[2] += ecross[l][p] * ecross[l][p];
    }

    if(iv != IntVect{0, 0, 0})
    {
        bool plc = (std::abs(conds[0] / 2. - 1.) > tolerance 
                    || std::abs(conds[2] / 2. - 1.) > tolerance);
        bool ppc = (std::abs(conds[1]) > tolerance);
        if (plc || ppc)
        {
            Print() << "---------\nLocation: " << iv << "\n";
            Print() << "Dot products: " << conds[0] / 2. << ", " << conds[2] / 2. << "\n";
            Print() << "Cross product: " << conds[1] << "\n";
            Print() << "Base tensor components: \n";
            for (int l=0; l<3; l++) for (int p=0; p<3; p++)
            {
                Print() << l << ", " << p << ": ";
                Print() << eplus[l][p] << ", " << ecross[l][p] << "\n";
            }
            Error("RandomField::Test_polarisation_tensor_orthonorm: polarisation tensors are not orthonormal here");
        }
    }
}

inline Real RandomField::calculate_total_power(const cMultiFab& fk, const int m_N)
{
    ReduceOps<ReduceOpSum> reduce_op;
    ReduceData<Real> reduce_data(reduce_op);
    using ReduceTuple = typename decltype(reduce_data)::Type;

    // 2. Loop over the grids owned by this MPI rank
    for (MFIter mfi(fk); mfi.isValid(); ++mfi) 
    {
        const Box& bx = mfi.fabbox();
        auto const& arr = fk.const_array(mfi);

        // 3. Execute the parallel reduction on the CPU/GPU
        reduce_op.eval(bx, reduce_data,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
            {
                // Get the real and imaginary parts
                Real re = arr(i,j,k).real();
                Real im = arr(i,j,k).imag();

                Real pow = re * re + im * im;
                // Multiply by 2 in most of the bulk,
                // to account for Hermitian modes.
                if (i != 0 && i != m_N/2) { pow *= 2.; }
                
                return pow;
            });
    }

    // 4. Extract the aggregated sum for this specific MPI rank
    ReduceTuple hv = reduce_data.value();
    Real total_local_power = get<0>(hv);

    // 5. Perform an MPI reduction to sum across all compute nodes
    ParallelDescriptor::ReduceRealSum(total_local_power);

    return total_local_power;
}

inline void RandomField::Test_Parsevals_thm(const MultiFab &hx, const cMultiFab &hk, const int m_N)
{
    Real xsum = std::pow(hx.norm2(), 2.);
    xsum /= std::pow(m_N, 3.);

    Real ksum = calculate_total_power(hk, m_N);

    int p = std::round(std::log10((ksum + ksum) / 2.));
    Real tol = tolerance * std::pow(10., p);

    if (std::abs(xsum - ksum) > tol)
    {
        Print() << "Order of magnitude of sum: " << p << "\n";
        Print() << "Tolerance: " << tol << "\n";
        Print() << "Stdev (x): " << xsum << "\n";
        Print() << "Integrated power (k): " << ksum << "\n";
        Print() << "Ratio: " << ksum / xsum << "\n";
        Print() << "Difference: " << std::abs(ksum - xsum) << "\n";
        Error("RandomField::Test_Parsevals_thm, Parseval's theorem fails here.");
    }
}

/****
    Device-callable computation functions
****/

// Returns analytic power spectrum in modulus/argument form
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
GpuComplex<Real> RandomField::calculate_mode_function(const Real km, const int spec_indx,
                                                       const Real H0)
{
    if(km < Real(1.e-23)) { return GpuComplex<Real>{0., 0.}; }

    Real ms_mag = 0.;
    Real ms_arg = 0.;

    Real kpr = km / H0;
    if (spec_indx == 0) // Position mode function
    {
        ms_mag = sqrt((1.0/km + H0*H0/pow(km, 3.))/2.);
        ms_arg = atan2((cos(kpr) + kpr*sin(kpr)), (kpr*cos(kpr) - sin(kpr)));
    }
    else if (spec_indx == 1) // Velocity mode function
    {
        ms_mag = sqrt(km/2.);
        ms_arg = -atan2(cos(kpr), sin(kpr));
    }
    else
    {
        AMREX_ASSERT_WITH_MESSAGE(false,
            "RandomField::calculate_mode_function Value of spec_type not allowed.");
    }

    return GpuComplex<Real>(ms_mag * cos(ms_arg), ms_mag * sin(ms_arg));
}

// Device-compatible STOIIC table lookup using pre-flattened device arrays
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
GpuComplex<Real> RandomField::find_in_stoiic_device(const Real km,
                                                      const int field_indx,
                                                      const FieldType field_type,
                                                      const GpuParams& gp)
{
    if(km == Real(0.)) { return GpuComplex<Real>{0., 0.}; }

    int spec_index = -1;
    for(int idx = 0; idx < gp.init_k_size; idx++)
    {
        if(std::abs(km - gp.d_init_k[idx]) < Real(1e-10))
        {
            spec_index = idx;
            break;
        }
    }
    AMREX_ASSERT_WITH_MESSAGE(spec_index >= 0,
        "RandomField::find_in_stoiic_device k not found in STOIIC table.");

    const Real* ps = (field_type == FieldType::Tensor) ? gp.d_tensor_ps : gp.d_scalar_ps;
    return GpuComplex<Real>{ps[(2*field_indx    )*gp.ps_row_size + spec_index],
                            ps[(2*field_indx + 1)*gp.ps_row_size + spec_index]};
}

// Host-only STOIIC lookup (accesses m_params Vectors)
// inline GpuComplex<Real> RandomField::find_in_stoiic(const Real km, const int field_indx,
//                                                       const FieldType field_type)
// {
//     GpuComplex<Real> zero(0., 0.);
//     if(km == 0) { return zero; }

//     int spec_index = -1;
//     for(int idx = 0; idx < (int)m_params.init_k.size(); idx++)
//     {
//         if(std::abs(km - m_params.init_k[idx]) < 1e-10)
//         {
//             spec_index = idx;
//             break;
//         }
//         else if (idx == (int)m_params.init_k.size() - 1)
//         {
//             Print() << km << "\n";
//             Error("RandomField::find_in_stoiic, "
//                   "The above k was not found in the STOIIC file.");
//         }
//     }

//     if(field_type == FieldType::Tensor)
//     {
//         return GpuComplex<Real>{m_params.tensor_ps[2*field_indx][spec_index],
//                                 m_params.tensor_ps[2*field_indx+1][spec_index]};
//     }
//     else if(field_type == FieldType::Scalar)
//     {
//         return GpuComplex<Real>{m_params.scalar_ps[2*field_indx][spec_index],
//                                 m_params.scalar_ps[2*field_indx+1][spec_index]};
//     }
//     else
//     {
//         Error("RandomField::find_in_stoiic field cannot be found.");
//         return GpuComplex<Real>{0., 0.};
//     }
// }

// Turns analytic PS into GRF and applies window function if requested
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
GpuComplex<Real> RandomField::calculate_random_field(const IntVect iv,
                                                      const int field_index,
                                                      const Real rand_amp,
                                                      const Real rand_phase,
                                                      const FieldType field_type,
                                                      const GpuParams& gp)
{
    GpuComplex<Real> value(0., 0.);

    Real kmag = get_kmag(iv, gp.Ni, gp.L);

    if(gp.read_from_stoiic)
    {
        value = find_in_stoiic_device(kmag, field_index, field_type, gp);
    }
    else
    {
        value = calculate_mode_function(kmag, field_index, gp.H0);
    }

    if(gp.use_rand == 1)
    {
        Real rand_mod = sqrt(-2. * log(rand_amp)); // Rayleigh distribution about |h|
        Real rand_arg = 2. * M_PI * rand_phase;    // Uniform random phase

        value *= rand_mod;

        Real new_real = value.real() * cos(rand_arg) - value.imag() * sin(rand_arg);
        Real new_imag = value.real() * sin(rand_arg) + value.imag() * cos(rand_arg);
        value = GpuComplex<Real>(new_real, new_imag);
    }

    if(gp.use_window == 1)
    {
        Real ks = (gp.N_coarse != 0 ?
                   std::sqrt(Real(3.)) * gp.N_coarse * M_PI / gp.L / 5. / 2. :
                   std::sqrt(Real(3.)) * gp.Ni       * M_PI / gp.L / 5. / 2.);
        Real Dt = gp.L / gp.Delta;
        value *= 0.5 * (1.0 - tanh(Dt * (kmag - ks)));
    }

    return value;
}

// Calculates basis vectors required for polarisation tensors
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
GpuArray<Real, 3> RandomField::calculate_basis_vector(const IntVect iv,
                                                      const int which_vector,
                                                      const int m_N,
                                                      const Real alpha)
{
    const Real i = static_cast<Real>(iv[0]);
    const Real j = static_cast<Real>(invert_index_with_sign(iv[1], m_N));
    const Real k = static_cast<Real>(invert_index_with_sign(iv[2], m_N));

    GpuArray<Real, 3> mhat{0., 0., 0.};
    GpuArray<Real, 3> nhat{0., 0., 0.};

    if (iv == IntVect{0, 0, 0}) { /* zero mode: return zeros */ }

    else if (i > 0.)
    {
        if (k == 0. && j == 0.)
        {
            mhat = {0., 1., 0.};
            nhat = {0., 0., 1.};
        }
        else
        {
            Real rnorm = sqrt((i*i + j*j) * (i*i + j*j + k*k));
            mhat = {j/sqrt(i*i + j*j), -i/sqrt(i*i + j*j), 0.};
            nhat = {(i*k) / rnorm, (j*k) / rnorm, -(i*i + j*j) / rnorm};
        }
    }

    else if (std::abs(j) > 0)
    {
        if(k == 0.)
        {
            mhat = {0., 0., 1.};
            nhat = {1., 0., 0.};
        }
        else
        {
            mhat = {-1., 0., 0.};
            nhat = {0., -k / sqrt(j*j + k*k), j / sqrt(j*j + k*k)};
        }
    }

    else if (std::abs(k) > 0)
    {
        mhat = {1., 0., 0.};
        nhat = {0., 1., 0.};
    }
    else
    {
        AMREX_ASSERT_WITH_MESSAGE(false,
            "RandomField::calculate_basis_vector Part of Fourier grid not covered.");
    }

    if (alpha != Real(0.))
    {
        Real a = alpha * M_PI / 180.;
        GpuArray<Real, 3> mp{0., 0., 0.};
        GpuArray<Real, 3> np{0., 0., 0.};
        for(int l = 0; l < 3; l++)
        {
            mp[l] =  cos(a) * mhat[l] + sin(a) * nhat[l];
            np[l] = -sin(a) * mhat[l] + cos(a) * nhat[l];
        }

        AMREX_ASSERT_WITH_MESSAGE(which_vector == 0 || which_vector == 1,
            "RandomField::calculate_basis_vector Incompatible vector type.");
        return (which_vector == 0) ? mp : np;
    }
    else
    {
        AMREX_ASSERT_WITH_MESSAGE(which_vector == 0 || which_vector == 1,
            "RandomField::calculate_basis_vector Incompatible vector type.");
        return (which_vector == 0) ? mhat : nhat;
    }
}

/****
    Nyquist conditions
****/

inline void RandomField::apply_nyquist_conditions(cMultiFab &field, const int m_N)
{
    const int nc = field.nComp();
    for (MFIter mfi(field); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.fabbox();
        Array4<GpuComplex<Real>> const& field_ptr = field.array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            if ((i == 0 || i == m_N/2) && (j == 0 || j == m_N/2) && (k == 0 || k == m_N/2))
            {
                for(int comp = 0; comp < nc; comp++)
                {
                    GpuComplex<Real> temp(field_ptr(i, j, k, comp).real(), 0.);
                    field_ptr(i, j, k, comp) = temp;
                }
            }
            else if (i == 0 || i == m_N/2)
            {
                if((k > m_N/2 && j == m_N/2) || (k == 0 && j > m_N/2) ||
                   (k > m_N/2 && j == 0)     || (k == m_N/2 && j > m_N/2))
                {
                    for(int comp = 0; comp < nc; comp++)
                    {
                        GpuComplex<Real> temp(
                             field_ptr(i, invert_index(j, m_N), invert_index(k, m_N), comp).real(),
                            -field_ptr(i, invert_index(j, m_N), invert_index(k, m_N), comp).imag());
                        field_ptr(i, j, k, comp) = temp;
                    }
                }
                else if(j > m_N/2)
                {
                    for(int comp = 0; comp < nc; comp++)
                    {
                        GpuComplex<Real> temp(
                             field_ptr(i, invert_index(j, m_N), flip_index(k, m_N), comp).real(),
                            -field_ptr(i, invert_index(j, m_N), flip_index(k, m_N), comp).imag());
                        field_ptr(i, j, k, comp) = temp;
                    }
                }
            }
        });
    }
}

/****
    Main initialisation routine
****/

inline void RandomField::init(amrex::MultiFab &state)
{
    BL_PROFILE("RandomField::init");

    // Derive MultiFab ingredients from state (configuration space)
    BoxArray sba = state.boxArray();
    DistributionMapping sdm = state.DistributionMap();

    int Ni = N;
    int dN = 1;
    if(m_params.N_fine != 0) { Ni = m_params.N_fine; dN = m_params.N_fine / N; }

    // Set up the problem domain in Fourier space
    // And impose that MPI ranks only slice along the i index (for Nyquist conditions)
    IntVect domain_low(0, 0, 0);
    BoxArray xba = (m_params.N_fine != 0 ? sba.refine(dN) : sba);

    IntVect k_domain_high(Ni/2, Ni-1, Ni-1);
    Box k_domain(domain_low, k_domain_high);
    Array< bool, AMREX_SPACEDIM > const &slicing{true, false, false};
    BoxArray kba = decompose(k_domain, ParallelContext::NProcsAll(), slicing);
    DistributionMapping kdm(kba);

    // Set up the MFs to store the in/out data sets
    cMultiFab hs_k(kba, kdm, 2, 0);
    cMultiFab As_k(kba, kdm, 2, 0);
    cMultiFab hij_k(kba, kdm, 6, 0);
    cMultiFab Aij_k(kba, kdm, 6, 0);
    MultiFab hij_x(xba, sdm, 6, 0);
    MultiFab Aij_x(xba, sdm, 6, 0);

    cMultiFab scalar_fields_k(kba, kdm, 4, 0);
    MultiFab scalar_fields_x(xba, sdm, 4, 0);

    hs_k.setVal(0.0);
    As_k.setVal(0.0);
    hij_k.setVal(0.0);
    Aij_k.setVal(0.0);
    hij_x.setVal(0.0);
    Aij_x.setVal(0.0);
    scalar_fields_k.setVal(0.0);
    scalar_fields_x.setVal(0.0);

    // Construct the Fourier transform
    IntVect x_domain_high(Ni-1, Ni-1, Ni-1);
    Box x_domain(domain_low, x_domain_high);
    FFT::R2C<Real> tensor_fft(x_domain, FFT::Info().setBatchSize(hij_k.nComp()));
    FFT::R2C<Real> scalar_fft(x_domain, FFT::Info().setBatchSize(scalar_fields_k.nComp()));

    // Build device parameter bundle — one value-capture replaces all this->X accesses
    GpuParams gp = make_gpu_params();
    if(m_params.N_fine != 0) { gp.Ni = Ni; }
    
    // Flatten STOIIC tables to device-accessible memory if needed
    Gpu::DeviceVector<Real> dv_init_k, dv_tensor_ps, dv_scalar_ps;
    if (m_params.read_from_stoiic)
    {
        dv_init_k = Gpu::DeviceVector<Real>(m_params.init_k.begin(), m_params.init_k.end());

        const int ps_row = static_cast<int>(m_params.init_k.size());
        Gpu::HostVector<Real> h_tensor(m_params.tensor_ps.size() * ps_row, 0.);
        for (int r = 0; r < (int)m_params.tensor_ps.size(); r++)
            for (int c = 0; c < ps_row; c++)
                h_tensor[r * ps_row + c] = m_params.tensor_ps[r][c];
        dv_tensor_ps = Gpu::DeviceVector<Real>(h_tensor.begin(), h_tensor.end());

        Gpu::HostVector<Real> h_scalar(m_params.scalar_ps.size() * ps_row, 0.);
        for (int r = 0; r < (int)m_params.scalar_ps.size(); r++)
            for (int c = 0; c < ps_row; c++)
                h_scalar[r * ps_row + c] = m_params.scalar_ps[r][c];
        dv_scalar_ps = Gpu::DeviceVector<Real>(h_scalar.begin(), h_scalar.end());

        gp.d_init_k    = dv_init_k.data();
        gp.d_tensor_ps = dv_tensor_ps.data();
        gp.d_scalar_ps = dv_scalar_ps.data();
        gp.init_k_size = static_cast<int>(dv_init_k.size());
        gp.ps_row_size = static_cast<int>(m_params.init_k.size());
    }

    Print() << "RandomField::init, Starting initial condition generation...\n";
    for (MFIter mfi(hs_k); mfi.isValid(); ++mfi) 
    {
        // Define the domain on this MPI rank
        const Box& bx = mfi.fabbox();

        // Make a pointer to the mode functions at this MPI box
        Array4<GpuComplex<Real>> const& hs_ptr = hs_k.array(mfi);
        Array4<GpuComplex<Real>> const& hij_ptr = hij_k.array(mfi);

        Array4<GpuComplex<Real>> const& As_ptr = As_k.array(mfi);
        Array4<GpuComplex<Real>> const& Aij_ptr = Aij_k.array(mfi);

        Array4<GpuComplex<Real>> const& scalar_fields_ptr = scalar_fields_k.array(mfi);

        // Loop to create mode functions, then hij(k) and Aij(k)
        amrex::ParallelFor(bx, 
            [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            IntVect iv = {i, j, k};
            amrex::InitRandom(gp.random_seed * (645950 * uint64_t(iv[0])
                                              + 520666 * uint64_t(invert_index_with_sign(iv[1], gp.Ni))
                                              + 767051 * uint64_t(invert_index_with_sign(iv[2], gp.Ni))
                                             ));

            if(gp.scalar_init)
            {
                Real draw1 = amrex::Random();
                Real draw2 = amrex::Random();

                for(int f = 0; f < 4; f++)
                {
                    scalar_fields_ptr(i, j, k, f) = calculate_random_field(
                        iv, f, draw1, draw2, FieldType::Scalar, gp);
                    if (gp.A != Real(1.0))
                    {
                        scalar_fields_ptr(i, j, k, f) *= gp.A;
                    }
                }
            }

            if(gp.tensor_init)
            {
                for(int p = 0; p < 2; p++)
                {
                    Real draw1 = amrex::Random();
                    Real draw2 = amrex::Random();

                    hs_ptr(i, j, k, p) = calculate_random_field(
                        iv, 0, draw1, draw2, FieldType::Tensor, gp);
                    As_ptr(i, j, k, p) = calculate_random_field(
                        iv, 1, draw1, draw2, FieldType::Tensor, gp);
                }

                GpuArray<Real, 3> mhat = calculate_basis_vector(iv, 0, gp.Ni, gp.alpha);
                GpuArray<Real, 3> nhat = calculate_basis_vector(iv, 1, gp.Ni, gp.alpha);

                Tensor<2, Real> eplus, ecross;
                for (int l = 0; l < 3; l++) for (int p = 0; p < 3; p++)
                {
                    eplus[l][p]  = mhat[l]*mhat[p] - nhat[l]*nhat[p];
                    ecross[l][p] = mhat[l]*nhat[p]  + nhat[l]*mhat[p];

                    hij_ptr(i, j, k, gp.lut[l][p]) =
                        hs_ptr(i, j, k, 0) * eplus[l][p] + hs_ptr(i, j, k, 1) * ecross[l][p];
                    Aij_ptr(i, j, k, gp.lut[l][p]) =
                        As_ptr(i, j, k, 0) * eplus[l][p] + As_ptr(i, j, k, 1) * ecross[l][p];
                }

                if (m_params.alpha != 0) { Test_polarisation_tensor_orthonorm(iv, eplus, ecross); }
            }
        });
    }

    // Apply the DC and Nyquist symmetry conditions
    apply_nyquist_conditions(hij_k, Ni);
    apply_nyquist_conditions(Aij_k, Ni);
    apply_nyquist_conditions(scalar_fields_k, Ni);

    // Do the Fourier transform
    tensor_fft.backward(hij_k, hij_x);
    tensor_fft.backward(Aij_k, Aij_x);
    scalar_fft.backward(scalar_fields_k, scalar_fields_x);

    // Apply normalisation into physical units
    // TODO: how does downsampling change this norm?
    hij_x.mult(norm);
    Aij_x.mult(norm);
    scalar_fields_x.mult(norm);

    Print() << "RandomField::init, Precision lost in phi is ";
    Print() << find_precision_loss(scalar_fields_x, 0, phi0) << "\n";
    Print() << "RandomField::init, Precision lost in chi is ";
    Print() << find_precision_loss(scalar_fields_x, 2, 1.0) << "\n";

    // Test that the resuling tensor perturbation field is trace-free
    Test_is_trace_free(hij_x, gp);

    // Convert to BSSN variables using the BSSN-CPT dictionary
    for (int l = 0; l < 3; l++) { hij_x.plus(1., lut[l][l], 1); }
    Aij_x.mult(-0.5);

    // Put these initial conditions into the state MF
    const int local_N = N;
    for (MFIter mfi(state); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.fabbox();
        Array4<Real> const& state_ptr  = state.array(mfi);
        Array4<Real> const& hij_ptr    = hij_x.array(mfi);
        Array4<Real> const& Aij_ptr    = Aij_x.array(mfi);
        Array4<Real> const& scalar_ptr = scalar_fields_x.array(mfi);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            const IntVect iv_ds{i, j, k};
            const IntVect iv{i * dN, j * dN, k * dN};

            if (iv_ds.min() >= 0 && iv_ds.max() < local_N)
            {
                if(gp.scalar_init)
                {
                    state_ptr(iv_ds, c_phi) += scalar_ptr(iv, 0);
                    state_ptr(iv_ds, c_Pi)  += scalar_ptr(iv, 1);
                    state_ptr(iv_ds, c_chi) += scalar_ptr(iv, 2);
                    state_ptr(iv_ds, c_K)   += scalar_ptr(iv, 3);
                }

                if(gp.tensor_init)
                {
                    state_ptr(iv_ds, c_h11) = hij_ptr(iv, gp.lut[0][0]);
                    state_ptr(iv_ds, c_h12) = hij_ptr(iv, gp.lut[0][1]);
                    state_ptr(iv_ds, c_h13) = hij_ptr(iv, gp.lut[0][2]);
                    state_ptr(iv_ds, c_h22) = hij_ptr(iv, gp.lut[1][1]);
                    state_ptr(iv_ds, c_h23) = hij_ptr(iv, gp.lut[1][2]);
                    state_ptr(iv_ds, c_h33) = hij_ptr(iv, gp.lut[2][2]);

                    state_ptr(iv_ds, c_A11) = Aij_ptr(iv, gp.lut[0][0]);
                    state_ptr(iv_ds, c_A12) = Aij_ptr(iv, gp.lut[0][1]);
                    state_ptr(iv_ds, c_A13) = Aij_ptr(iv, gp.lut[0][2]);
                    state_ptr(iv_ds, c_A22) = Aij_ptr(iv, gp.lut[1][1]);
                    state_ptr(iv_ds, c_A23) = Aij_ptr(iv, gp.lut[1][2]);
                    state_ptr(iv_ds, c_A33) = Aij_ptr(iv, gp.lut[2][2]);
                }
            }
        });
    }
}

/****
    Extraction routines
****/

// Calculates and prints the power spectrum.
// Uses Gpu::DeviceVector for ps_map and kcount to allow device atomic adds.
inline void RandomField::print_power_spectrum(cMultiFab &field_array, SmallDataIO &power_spec_file, const int component)
{
    const int local_N   = N;
    const Real kiso_max = std::sqrt(3.) * local_N * M_PI / m_params.L;
    const Real dkiso    = sqrt(3.) * 2. * M_PI / m_params.L;
    const Real tol      = 1.e-12;

    // check the stepping along the diagonal is consistent
    if (kiso_max/dkiso - local_N/2 > tol)
    {
        Error("RandomField::print_power_spectrum Isotropic k axis is too large.");
    }

    Vector<Real> kiso(local_N/2 + 1, 0.);
    for (int s = 0; s <= local_N/2; s++) { kiso[s] = s * dkiso; }

    // Device-accessible accumulation arrays (Issue 6 fix: no reference capture)
    Gpu::DeviceVector<Real> d_ps_map(local_N/2 + 1, 0.);
    Gpu::DeviceVector<int>  d_kcount(local_N/2 + 1, 0);
    Real* d_ps_ptr     = d_ps_map.data();
    int*  d_kcount_ptr = d_kcount.data();

    // kiso array on device
    Gpu::DeviceVector<Real> d_kiso(kiso.begin(), kiso.end());
    const Real* d_kiso_ptr = d_kiso.data();

    GpuParams gp = make_gpu_params();

    for (MFIter mfi(field_array); mfi.isValid(); ++mfi)
    {
        Array4<GpuComplex<Real>> const& field_ptr = field_array.array(mfi);
        const Box& bx = mfi.fabbox();

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            IntVect iv{i, j, k};
            Real kmag = get_kmag(iv, local_N, gp.L);

            AMREX_ASSERT_WITH_MESSAGE(!(kmag - kiso_max > tol),
                "RandomField::print_power_spectrum Found magnitude larger than (N/2,N/2,N/2).");

            for (int s = 1; s <= local_N/2; s++)
            {
                AMREX_ASSERT_WITH_MESSAGE(!(kmag < d_kiso_ptr[0]),
                    "RandomField::print_power_spectrum kmag below the kiso domain.");
                AMREX_ASSERT_WITH_MESSAGE(!(kmag - d_kiso_ptr[local_N/2] > tol),
                    "RandomField::print_power_spectrum kmag above the kiso domain.");

                if (kmag < d_kiso_ptr[s] && kmag >= d_kiso_ptr[s-1])
                {
                    Real power = (pow(field_ptr(i, j, k, component).real(), 2.)
                                + pow(field_ptr(i, j, k, component).imag(), 2.));
                    if (i != 0 && i != local_N/2) { power *= 2.; }

                    Gpu::Atomic::Add(d_kcount_ptr + (s-1), 1);
                    Gpu::Atomic::Add(d_ps_ptr + (s-1), power);
                    break;
                }
                else if (kmag == d_kiso_ptr[local_N/2])
                {
                    Real power = (pow(field_ptr(i, j, k, component).real(), 2.)
                                + pow(field_ptr(i, j, k, component).imag(), 2.));
                    if (i != 0 && i != local_N/2) { power *= 2.; }

                    Gpu::Atomic::Add(d_kcount_ptr + local_N/2, 1);
                    Gpu::Atomic::Add(d_ps_ptr + local_N/2, power);
                    break;
                }
                else if(s > local_N/2)
                {
                    AMREX_ASSERT_WITH_MESSAGE(false,
                        "RandomField::print_power_spectrum Part of the spectrum isn't captured.");
                }
            }
        });
    }

    // Copy results back to host
    Vector<Real> ps_map(local_N/2 + 1, 0.);
    Vector<int>  kcount(local_N/2 + 1, 0);
    Gpu::copy(Gpu::deviceToHost, d_ps_map.begin(), d_ps_map.end(), ps_map.begin());
    Gpu::copy(Gpu::deviceToHost, d_kcount.begin(), d_kcount.end(), kcount.begin());

    ParallelAllReduce::Sum(kcount.data(),  static_cast<int>(kcount.size()),  ParallelContext::CommunicatorSub());
    ParallelAllReduce::Sum(ps_map.data(),  static_cast<int>(ps_map.size()),  ParallelContext::CommunicatorSub());

    // Print the power spectrum to a new file in data/
#pragma omp single
    for(int s = 0; s < N/2; s++)
    {
        power_spec_file.write_data_line({(kiso[s]+kiso[s+1])/2., (Real)ps_map[s]/kcount[s]});
    }
}

// Finds statistical moment x of given MultiFab using ReduceOps (Issue 6 fix)
inline Real RandomField::find_field_moment_x(MultiFab &field, const Vector<Real> mean,
                                              const int moment, const int component)
{
    const Real vol = std::pow(N, 3.);
    const Real mean_comp = mean[component];

    ReduceOps<ReduceOpSum> reduce_op;
    ReduceData<Real> reduce_data(reduce_op);
    using ReduceTuple = typename decltype(reduce_data)::Type;

    for (MFIter mfi(field); mfi.isValid(); ++mfi)
    {
        Array4<Real> const& field_ptr = field.array(mfi);
        const Box& bx = mfi.fabbox();

        reduce_op.eval(bx, reduce_data,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
            {
                return { pow(field_ptr(i, j, k, component) - mean_comp, moment) };
            });
    }

    Real sum = get<0>(reduce_data.value());
    ParallelAllReduce::Sum(sum, ParallelContext::CommunicatorSub());

    // Normalise and return moment x
    if (sum == 0) { return 0; }
    else if(moment == 2) { return sqrt(sum/vol); }
    else                 { return sum/vol; }
}

// Calculates and prints requested moments (any between 1 and 4)
inline Vector<Real> RandomField::print_moment(MultiFab &field, const Vector<std::string> names,
                                              const Vector<int> &moment_orders, SmallDataIO &file,
                                              const int is_first_step)
{
    // Trap instance where the user requests too large a moment
    for(const auto moment : moment_orders)
    {
        if(moment > 4)
        {
            Error("RandomField::print_moment Chosen moment order has not been implemented");
        }
    }

    // Allocate arrays to store each moment
    const int nc = field.nComp();
    const Real vol = std::pow(N, 3.);
    Vector<Real> means(nc, 0.);
    Vector<Real> stdev(nc, 0.);
    Vector<Real> skew(nc, 0.);
    Vector<Real> kurt(nc, 0.);

    // Find iterators, which determine which moments are requested and their ordering
    Vector<int>::const_iterator start = moment_orders.begin();
    Vector<int>::const_iterator mean_itr = std::find(moment_orders.begin(), moment_orders.end(), 1);
    Vector<int>::const_iterator stdev_itr= std::find(moment_orders.begin(), moment_orders.end(), 2);
    Vector<int>::const_iterator skew_itr = std::find(moment_orders.begin(), moment_orders.end(), 3);
    Vector<int>::const_iterator kurt_itr = std::find(moment_orders.begin(), moment_orders.end(), 4);

    // Allocate vectors to store header line and data lines
    Vector<Real> data_to_print(nc * moment_orders.size(), 0.);
    Vector<std::string> headers(nc * moment_orders.size(), "");

    for (int comp = 0; comp < nc; comp++)
    {
        means[comp] = field.sum(comp)/vol;
        if(mean_itr != moment_orders.end())
        {
            assign_statistics_data(headers, names[comp]+"-mean", data_to_print, means, comp, nc,
                                    mean_itr, start, is_first_step);
        }

        if(moment_orders.back() != 1)
        {
            stdev[comp] = find_field_moment_x(field, means, 2, comp);
            if(stdev_itr != moment_orders.end())
            {
                assign_statistics_data(headers, names[comp]+"-stdev", data_to_print, stdev, comp, nc,
                                        stdev_itr, start, is_first_step);
            }

            if(moment_orders.back() != 2)
            {
                skew[comp] = find_field_moment_x(field, means, 3, comp);
                skew[comp] /= std::pow(stdev[comp], 3.);

                if(skew_itr != moment_orders.end())
                {
                    assign_statistics_data(headers, names[comp]+"-skew", data_to_print, skew, comp, nc,
                                            skew_itr, start, is_first_step);
                }

                if(moment_orders.back() != 3)
                {
                    kurt[comp] = find_field_moment_x(field, means, 4, comp);
                    kurt[comp] /= std::pow(stdev[comp], 4.);

                    assign_statistics_data(headers, names[comp]+"-kurt", data_to_print, kurt, comp, nc,
                                            kurt_itr, start, is_first_step);
                }
            }
        }
    }

    // Write header and data lines
#pragma omp single
    if(is_first_step) { file.write_header_line(headers); }

#pragma omp single
    file.write_time_data_line(data_to_print);

    return stdev;
}

inline void RandomField::derive(const MultiFab &state, MultiFab &out, int dcomp)
{
    BL_PROFILE("RandomField::derive");

    // Extract MultiFab ingredients from state
    BoxArray sba = state.boxArray();
    DistributionMapping sdm = state.DistributionMap();
    MultiFab gij_x(sba, sdm, 6, 0);

    // 0: scalar field
    // 1: conformal factor
    MultiFab scalars_x(sba, sdm, 2, 0);

    // Copy the spatial metric from the state
    Copy(gij_x, state, c_h11, lut[0][0], 1, 0);
    Copy(gij_x, state, c_h12, lut[0][1], 1, 0);
    Copy(gij_x, state, c_h13, lut[0][2], 1, 0);
    Copy(gij_x, state, c_h22, lut[1][1], 1, 0);
    Copy(gij_x, state, c_h23, lut[1][2], 1, 0);
    Copy(gij_x, state, c_h33, lut[2][2], 1, 0);

    const int m_c_phi = 0;
    const int m_c_chi = 1;
    Copy(scalars_x, state, c_phi, m_c_phi, 1, 0);
    Copy(scalars_x, state, c_chi, m_c_chi, 1, 0);

    // Find background quantities needed to extract \cal R
    const int vol = std::pow(m_params.N_readin, 3);
    const Real K_bar = state.sum(c_K)/vol;
    const Real alpha_bar = state.sum(c_lapse)/vol;
    const Real Pi_bar = state.sum(c_Pi)/vol;
    const Real phi_bar = state.sum(c_phi)/vol;
    const Real chi_bar = state.sum(c_chi)/vol;

    // Remove background from scalar field
    scalars_x.plus(-phi_bar, m_c_phi, 1);
    scalars_x.plus(-chi_bar, m_c_chi, 1);

    for (int l = 0; l < 3; l++) { gij_x.plus(-1., lut[l][l], 1); }

    // Set up the problem domain in Fourier space
    // And impose that MPI ranks only slice along the i index (for Nyquist conditions)
    IntVect domain_low(0, 0, 0);
    IntVect k_domain_high(N/2, N-1, N-1);
    Box k_domain(domain_low, k_domain_high);
    Array< bool, AMREX_SPACEDIM > const &slicing{true, false, false};
    BoxArray kba = decompose(k_domain, ParallelContext::NProcsAll(), slicing);
    DistributionMapping kdm(kba);

    // Set up the arrays to store the Fourier data sets
    cMultiFab hs_k(kba, kdm, 2, 0);
    cMultiFab gij_k(kba, kdm, 6, 0);
    cMultiFab scalars_k(kba, kdm, 2, 0);
    cMultiFab R_k(kba, kdm, 1, 0);

    hs_k.setVal(0.0);
    gij_k.setVal(0.0);
    scalars_k.setVal(0.0);
    R_k.setVal(0.0);

    // Set up the FFT
    IntVect x_domain_high(N-1, N-1, N-1);
    Box x_domain(domain_low, x_domain_high);
    FFT::R2C<Real> tensor_fft(x_domain, FFT::Info().setBatchSize(gij_k.nComp()));
    FFT::R2C<Real> scalar_fft(x_domain, FFT::Info().setBatchSize(scalars_k.nComp()));

    // Perform the fft
    tensor_fft.forward(gij_x, gij_k);
    scalar_fft.forward(scalars_x, scalars_k);

    // Normalise the fft (fftw style)
    for(int comp = 0; comp < 6; comp++) { gij_k.mult(1./norm/pow(N, 3.), comp, 1); }
    for(int comp = 0; comp < 2; comp++) { scalars_k.mult(1./norm/pow(N, 3.), comp, 1); }

    // Loop to extract the Fourier-space mode functions
    for (MFIter mfi(gij_k); mfi.isValid(); ++mfi) 
    {
        const Box& bx = mfi.fabbox();
        Array4<GpuComplex<Real>> const& hs_ptr      = hs_k.array(mfi);
        Array4<GpuComplex<Real>> const& hij_ptr      = gij_k.array(mfi);
        Array4<GpuComplex<Real>> const& scalars_ptr  = scalars_k.array(mfi);
        Array4<GpuComplex<Real>> const& R_k_ptr      = R_k.array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            IntVect iv{i, j, k};

            if (iv != IntVect{0, 0, 0})
            {
                GpuArray<Real, 3> mhat = calculate_basis_vector(iv, 0, gp.N, gp.alpha);
                GpuArray<Real, 3> nhat = calculate_basis_vector(iv, 1, gp.N, gp.alpha);
                Tensor<2, Real> eplus, ecross;

                for (int l = 0; l < 3; l++) for (int p = 0; p < 3; p++)
                {
                    eplus[l][p]  = mhat[l]*mhat[p] - nhat[l]*nhat[p];
                    ecross[l][p] = mhat[l]*nhat[p]  + nhat[l]*mhat[p];

                    hs_ptr(i, j, k, 0) += (hij_ptr(i, j, k, gp.lut[l][p]) * eplus[l][p]) / 2.;
                    hs_ptr(i, j, k, 1) += (hij_ptr(i, j, k, gp.lut[l][p]) * ecross[l][p]) / 2.;
                }

                if (m_params.alpha != 0) { Test_polarisation_tensor_orthonorm(iv, eplus, ecross); }

                // Calculate the TT and scalar-(vector) components of the 
                // metric, by reconstructing hij and subtracting it from \tilde{gamma}_ij
                Tensor<2, GpuComplex<Real>> hij, hSV;
                for (int l = 0; l < 3; l++) for (int p = 0; p < 3; p++)
                {
                    hij[l][p] = hs_ptr(i, j, k, 0) * eplus[l][p] + hs_ptr(i, j, k, 1) * ecross[l][p];
                    hSV[l][p] = hij_ptr(i, j, k, gp.lut[l][p]) - hij[l][p];
                }

                if(gp.scalar_init)
                {
                    GpuArray<Real, 3> iv_k = {
                        static_cast<Real>(iv[0]),
                        static_cast<Real>(invert_index_with_sign(iv[1], gp.N)),
                        static_cast<Real>(invert_index_with_sign(iv[2], gp.N))
                    };
                    for(int d = 0; d < 3; d++) { iv_k[d] *= 2. * M_PI / gp.L; }

                    Real kmag = get_kmag(iv, gp.N, gp.L);
                    GpuComplex<Real> Phi = 0;

                    if(kmag == Real(0.))
                    {
                        R_k_ptr(i, j, k) = GpuComplex<Real>{0., 0.};
                    }
                    else
                    {
                        for(int l = 0; l < 3; l++) for(int p = 0; p < 3; p++)
                        {
                            Phi += (iv_k[l] * iv_k[p] * hSV[l][p]) / std::pow(kmag, 2.);
                        }
                        Phi *= 1./4.;
                        Phi += 0.5 * scalars_ptr(i, j, k, m_c_chi);

                        R_k_ptr(i, j, k) = Phi - (K_bar/3.) * scalars_ptr(i, j, k, m_c_phi)
                                               / alpha_bar / Pi_bar;
                    }
                }
            }
        });
    }

    apply_nyquist_conditions(hs_k, N);
    apply_nyquist_conditions(R_k, N);

    // Make a multifab to store config space mode functions
    // Need to use out to make these ingredients??
    BoxArray xba = out.boxArray();//(x_domain); //
    DistributionMapping xdm = out.DistributionMap();//(xba); //
    MultiFab hs_x(xba, xdm, 2, 0);
    MultiFab R_x(xba, xdm, 1, 0);
    hs_x.setVal(0.0);
    R_x.setVal(0.0);

    // Fourier transform
    FFT::R2C<Real> mode_function_fft(x_domain, FFT::Info().setBatchSize(hs_k.nComp()));
    FFT::R2C<Real> R_fft(x_domain, FFT::Info().setBatchSize(R_x.nComp()));
    mode_function_fft.backward(hs_k, hs_x);
    R_fft.backward(R_k, R_x);

    Test_Parsevals_thm(hs_x, hs_k, N);
    Test_Parsevals_thm(R_x, R_k, N);

    hs_x.mult(norm);
    R_x.mult(norm);

    for (MFIter mfi(hs_x); mfi.isValid(); ++mfi) 
    {
        Array4<Real> const& out_ptr = out.array(mfi);
        Array4<Real> const& hx_ptr = hs_x.array(mfi);
        Array4<Real> const& Rx_ptr = R_x.array(mfi);

        const Box& bx = mfi.fabbox();
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            const IntVect iv{i, j, k};

            out_ptr(iv, dcomp) = hx_ptr(i, j, k, 0);
            out_ptr(iv, dcomp + 1) = hx_ptr(i, j, k, 1);
            out_ptr(iv, dcomp + 2) = Rx_ptr(i, j, k);
        });
    }
}

// Main extraction routine
inline void RandomField::extract(const MultiFab &state, const std::string data_path, const Real dt,  
                                 const Real cur_time, const int restart_time, const int first_step)
{
    BL_PROFILE("RandomField::extract");

    // Extract MultiFab ingredients from state
    BoxArray sba = state.boxArray();
    DistributionMapping sdm = state.DistributionMap();
    MultiFab gij_x(sba, sdm, 6, 0);

    // 0: scalar field
    // 1: conformal factor
    MultiFab scalars_x(sba, sdm, 2, 0);

    // Copy the spatial metric from the state
    Copy(gij_x, state, c_h11, lut[0][0], 1, 0);
    Copy(gij_x, state, c_h12, lut[0][1], 1, 0);
    Copy(gij_x, state, c_h13, lut[0][2], 1, 0);
    Copy(gij_x, state, c_h22, lut[1][1], 1, 0);
    Copy(gij_x, state, c_h23, lut[1][2], 1, 0);
    Copy(gij_x, state, c_h33, lut[2][2], 1, 0);

    const int m_c_phi = 0;
    const int m_c_chi = 1;
    Copy(scalars_x, state, c_phi, m_c_phi, 1, 0);
    Copy(scalars_x, state, c_chi, m_c_chi, 1, 0);

    // Find background quantities needed to extract \cal R
    const int vol = std::pow(N, 3);
    const Real K_bar = state.sum(c_K)/vol;
    const Real alpha_bar = state.sum(c_lapse)/vol;
    const Real Pi_bar = state.sum(c_Pi)/vol;
    const Real phi_bar = state.sum(c_phi)/static_cast<Real>(vol);
    const Real chi_bar = state.sum(c_chi)/static_cast<Real>(vol);

    // Remove background from scalar field
    scalars_x.plus(-phi_bar, m_c_phi, 1);
    scalars_x.plus(-chi_bar, m_c_chi, 1);

    // Undo the normalisation and BSSN-CPT conversion
    for (int l=0; l<3; l++) { gij_x.plus(-1., lut[l][l], 1); }

    // Set up the problem domain in Fourier space
    // And impose that MPI ranks only slice along the i index (for Nyquist conditions)
    IntVect domain_low(0, 0, 0);
    IntVect k_domain_high(N/2, N-1, N-1);
    Box k_domain(domain_low, k_domain_high);
    Array< bool, AMREX_SPACEDIM > const &slicing{true, false, false};
    BoxArray kba = decompose(k_domain, ParallelContext::NProcsAll(), slicing);
    DistributionMapping kdm(kba);

    // Set up the arrays to store the Fourier data sets
    cMultiFab hs_k(kba, kdm, 2, 0);
    cMultiFab gij_k(kba, kdm, 6, 0);
    cMultiFab scalars_k(kba, kdm, 2, 0);
    cMultiFab R_k(kba, kdm, 1, 0);

    hs_k.setVal(0.0);
    gij_k.setVal(0.0);
    scalars_k.setVal(0.0);
    R_k.setVal(0.0);

    // Set up the FFT
    IntVect x_domain_high(N-1, N-1, N-1);
    Box x_domain(domain_low, x_domain_high);
    FFT::R2C<Real> tensor_fft(x_domain, FFT::Info().setBatchSize(gij_k.nComp()));
    FFT::R2C<Real> scalar_fft(x_domain, FFT::Info().setBatchSize(scalars_k.nComp()));

    // Perform the fft
    tensor_fft.forward(gij_x, gij_k);
    scalar_fft.forward(scalars_x, scalars_k);

    // Normalise the fft (fftw style)
    for(int comp = 0; comp < 6; comp++) { gij_k.mult(1./norm/pow(N, 3.), comp, 1); }
    for(int comp = 0; comp < 2; comp++) { scalars_k.mult(1./norm/pow(N, 3.), comp, 1); }

    GpuParams gp = make_gpu_params();

    // Track max traces via ReduceOps (Issue 6 fix: no reference capture of host scalars)
    Real hij_tr_max = 0.;
    Real hSV_tr_max = 0.;

    // Loop to extract the Fourier-space mode functions
    for (MFIter mfi(gij_k); mfi.isValid(); ++mfi) 
    {
        const Box& bx = mfi.fabbox();
        Array4<GpuComplex<Real>> const& hs_ptr     = hs_k.array(mfi);
        Array4<GpuComplex<Real>> const& hij_ptr     = gij_k.array(mfi);
        Array4<GpuComplex<Real>> const& scalars_ptr = scalars_k.array(mfi);
        Array4<GpuComplex<Real>> const& R_k_ptr     = R_k.array(mfi);

        ReduceOps<ReduceOpMax, ReduceOpMax> reduce_op;
        ReduceData<Real, Real> reduce_data(reduce_op);
        using ReduceTuple = typename decltype(reduce_data)::Type;

        reduce_op.eval(bx, reduce_data,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
        {
            IntVect iv{i, j, k};

            GpuArray<Real, 3> mhat = calculate_basis_vector(iv, 0, gp.N, gp.alpha);
            GpuArray<Real, 3> nhat = calculate_basis_vector(iv, 1, gp.N, gp.alpha);
            Tensor<2, Real> eplus, ecross;

            // Find basis tensors and do the Fourier trick
            for (int l=0; l<3; l++) for (int p=0; p<3; p++)
            {
                eplus[l][p] = mhat[l]*mhat[p] - nhat[l]*nhat[p];
                ecross[l][p] = mhat[l]*nhat[p] + nhat[l]*mhat[p];

                hs_ptr(i, j, k, 0) += (hij_ptr(i, j, k, gp.lut[l][p]) * eplus[l][p]) / 2.;
                hs_ptr(i, j, k, 1) += (hij_ptr(i, j, k, gp.lut[l][p]) * ecross[l][p]) / 2.;
            }

            if (m_params.alpha != 0) { Test_polarisation_tensor_orthonorm(iv, eplus, ecross); }

            // Calculate the TT and scalar-(vector) components of the 
            // metric, by reconstructing hij and subtracting it from \tilde{gamma}_ij
            Tensor<2, GpuComplex<Real>> hij, hSV;
            GpuComplex<Real> hij_tr = 0.;
            GpuComplex<Real> hSV_tr = 0.;
            for (int l=0; l<3; l++) for (int p=0; p<3; p++)
            {
                hij[l][p] = hs_ptr(i, j, k, 0) * eplus[l][p] + hs_ptr(i, j, k, 1) * ecross[l][p];
                hSV[l][p] = hij_ptr(i, j, k, gp.lut[l][p]) - hij[l][p];
                if(l == p) { hij_tr += hij[l][p]; hSV_tr += hSV[l][p]; }
            }

            Real hij_tr_mag = sqrt(pow(hij_tr.real(), 2.) + pow(hij_tr.imag(), 2.));
            Real hSV_tr_mag = sqrt(pow(hSV_tr.real(), 2.) + pow(hSV_tr.imag(), 2.));

            AMREX_ASSERT_WITH_MESSAGE(!(std::abs(hij_tr_mag) > gp.tolerance),
                "RandomField::extract, hij trace magnitude too large in extraction");

            if(gp.scalar_init)
            {
                GpuArray<Real, 3> iv_k = {
                    static_cast<Real>(iv[0]),
                    static_cast<Real>(invert_index_with_sign(iv[1], gp.N)),
                    static_cast<Real>(invert_index_with_sign(iv[2], gp.N))
                };
                for(int d = 0; d < 3; d++) { iv_k[d] *= 2. * M_PI / gp.L; }

                Real kmag = get_kmag(iv, gp.N, gp.L);
                GpuComplex<Real> Phi = 0;

                if(kmag == Real(0.))
                {
                    R_k_ptr(i, j, k) = GpuComplex<Real>{0., 0.};
                }

                else
                {
                    // converstion from chi and gamma_ij -> Phi
                    for(int l=0; l<3; l++) for(int p=0; p<3; p++)
                    {
                        Phi += (iv_k[l] * iv_k[p] * hSV[l][p])/std::pow(kmag, 2.);
                    }
                    Phi *= 1./4.;
                    Phi += 0.5 * scalars_ptr(i, j, k, m_c_chi);
                    R_k_ptr(i, j, k) = Phi - (K_bar/3.) * scalars_ptr(i, j, k, m_c_phi)
                                           / alpha_bar / Pi_bar;
                }
            }

            return { hij_tr_mag, hSV_tr_mag };
        });

        ReduceTuple hv = reduce_data.value();
        hij_tr_max = amrex::max(hij_tr_max, amrex::get<0>(hv));
        hSV_tr_max = amrex::max(hSV_tr_max, amrex::get<1>(hv));
    }

    SmallDataIO trace_file(data_path+"tensor-traces", dt, cur_time,
                           restart_time, SmallDataIO::APPEND, first_step, ".dat");
    if(first_step) { trace_file.write_header_line({"hij trace max", "hSV trace max"}); }
    trace_file.write_time_data_line({hij_tr_max, hSV_tr_max});

    apply_nyquist_conditions(hs_k, N);
    apply_nyquist_conditions(R_k, N);

    // Find the binned PS for each mode function and print to data/
    if((m_params.calc_binned_power_spectrum) 
	    && (static_cast<int>(cur_time/dt) % m_params.plot_int == 0))
    {
	    Print() << "RandomField::extract, Time step at print: ";
        Print() << static_cast<int>(std::round(cur_time/dt)) << "\n";

        std::string spec_path = make_subdirectory(data_path, "spectra", first_step);
        Vector<std::string> filenames(2, "");

        for(int comp = 0; comp < hs_k.nComp(); comp++)
        {
            filenames[comp] = spec_path+"spectrum-comp-"+std::to_string(comp)+"-time-";
            SmallDataIO spectrum_file(filenames[comp], dt, cur_time, restart_time, SmallDataIO::NEW, first_step, ".dat");
            print_power_spectrum(hs_k, spectrum_file, comp);
        }

        std::string filename = spec_path+"spectrum-Rk-time-";
        SmallDataIO spectrum_file(filename, dt, cur_time, restart_time, SmallDataIO::NEW, first_step, ".dat");
        print_power_spectrum(R_k, spectrum_file, 0);
    }

    // Find mode functions in configuration space if requested,
    // and find the statistics (orders 1-4) of the polarisation fields 
    // and the R field. 
    // And print the tensor-to-scalar ratio if requested.
    if(m_params.calc_higher_order_statistics)
    {
        // Make a multifab to store config space mode functions
        BoxArray xba(x_domain);
        DistributionMapping xdm(xba);
        MultiFab hs_x(xba, xdm, 2, 0);
        MultiFab R_x(xba, xdm, 1, 0);
        hs_x.setVal(0.0);
        R_x.setVal(0.0);

        // Fourier transform
        FFT::R2C<Real> tensor_mode_function_fft(x_domain, FFT::Info().setBatchSize(hs_k.nComp()));
        FFT::R2C<Real> scalar_mode_function_fft(x_domain, FFT::Info().setBatchSize(R_k.nComp()));
        tensor_mode_function_fft.backward(hs_k, hs_x);
        scalar_mode_function_fft.backward(R_k, R_x);

        if (m_params.tensor_init) 
        { 
            for (int c=0; c<hs_x.nComp(); c++)
            {
                MultiFab hx_tmp(hs_x, amrex::make_alias, c, 1);
                cMultiFab hk_tmp(hs_k, amrex::make_alias, c, 1);
                Test_Parsevals_thm(hx_tmp, hk_tmp);
            } 
        }
        if (m_params.scalar_init) { Test_Parsevals_thm(R_x, R_k, N); }

        // Apply physical normalisation
        hs_x.mult(norm);
        R_x.mult(norm);

        int output_comps = hs_x.nComp() + R_x.nComp();
        MultiFab out_MF(hs_x.boxArray(), hs_x.DistributionMap(), output_comps, 0);
        Copy(out_MF, R_x, 0, 0, R_k.nComp(), 0);
        Copy(out_MF, hs_x, 0, R_k.nComp(), hs_x.nComp(), 0);

        // Print mode functions if requested
        if(m_params.print_mode_functions == 1)
        {
            std::string mf_path = make_subdirectory(data_path, "mode-functions", first_step);
            std::string mf_filename = mf_path + "mode-function-" + std::to_string(cur_time/dt);
            const int local_N = N;

            // File I/O must run on the host — use a CPU loop, not ParallelFor
            for (MFIter mfi(hs_x); mfi.isValid(); ++mfi)
            {
                Array4<Real> const& hx_ptr = hs_x.array(mfi);
                const Box& bx = mfi.fabbox();
                const auto lo = bx.smallEnd();
                const auto hi = bx.bigEnd();

                for (int kk = lo[2]; kk <= hi[2]; ++kk)
                for (int jj = lo[1]; jj <= hi[1]; ++jj)
                for (int ii = lo[0]; ii <= hi[0]; ++ii)
                {
                    AllPrintToFile(mf_filename) << ii + local_N*(jj + local_N*kk) << ",";
                    for(int c = 0; c < 2; c++)
                    {
                        AllPrintToFile(mf_filename).SetPrecision(14) << hx_ptr(ii, jj, kk, c) << ",";
                    }
                    AllPrintToFile(mf_filename) << "\n";
                }
            }
        }

        // Calculate and print field moments if requested
        Vector<Real> stdevs;
        if (m_params.calc_higher_order_statistics)
        {
            SmallDataIO stats_file(data_path+"field-statistics", dt, cur_time, 
                                    restart_time, SmallDataIO::APPEND, first_step, ".dat");

            if (!m_params.orders.empty())
            {
                Vector<std::string> names{"R","hplus","hcross"};
                stdevs = print_moment(out_MF, names, m_params.orders, stats_file, first_step);
                // Print() << stdevs[0] << "\n";
                // Error();
            }
        }
        
        SmallDataIO ts_file(data_path+"tensor-scalar-ratio", dt, cur_time, 
                                    restart_time, SmallDataIO::APPEND, first_step, ".dat");
#pragma omp single
        if(first_step) { ts_file.write_header_line({"T/S ratio (plus)", "T/S ratio (cross)"}); }

#pragma omp single
    	ts_file.write_time_data_line({stdevs[1] / stdevs[0], stdevs[2] / stdevs[0]}); 
    }
}

#endif /* RANDOMFIELD_IMPL_HPP_*/

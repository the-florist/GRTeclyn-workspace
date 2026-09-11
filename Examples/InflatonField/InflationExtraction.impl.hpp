/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#if !defined(INFLATIONEXTRACTION_HPP_)
#error "This file should only be included via InflationExtraction.hpp"
#endif

#ifndef INFLATIONEXTRACTION_IMPL_HPP_
#define INFLATIONEXTRACTION_IMPL_HPP_

#include <AMReX_ParallelReduce.H>
#include <AMReX_Reduce.H>

#include <algorithm>
#include <cmath>

/* Helper functions */

inline std::string
InflationExtraction::make_subdirectory(const std::string &subdirectory) const
{
    const std::string new_path = m_params.data_path + subdirectory + "/";
    // Idempotent, and cannot be restricted to the first step: a restarted run
    // writes here too
    FilesystemTools::ensure_directory_exists(new_path);
    return new_path;
}

inline void InflationExtraction::assign_statistics_data(
    amrex::Vector<std::string> &header_storage, const std::string &name,
    amrex::Vector<amrex::Real> &data_storage, const amrex::Real value,
    const int component, const int num_comps,
    const amrex::Vector<int>::const_iterator itr,
    const amrex::Vector<int>::const_iterator start, const bool is_first_step)
{
    const int loc = component + num_comps * static_cast<int>(itr - start);
    if (is_first_step)
    {
        header_storage[loc] = name;
    }
    data_storage[loc] = value;
}

// Volume averages of the evolved background quantities
inline void
InflationExtraction::extract_background_means(const amrex::MultiFab &state)
{
    // The unitless grid volume, matching the convention used for the field
    // moments
    const amrex::Real vol = std::pow(params().N, 3.);

    // MultiFab::sum reduces across ranks, so every rank sees the same average
    const amrex::Real phi_mean   = state.sum(c_phi) / vol;
    const amrex::Real pi_mean    = state.sum(c_Pi) / vol;
    const amrex::Real chi_mean   = state.sum(c_chi) / vol;
    const amrex::Real k_mean     = state.sum(c_K) / vol;
    const amrex::Real lapse_mean = state.sum(c_lapse) / vol;

    // The FLRW background these averages correspond to. Both are built from
    // the averaged variable rather than averaged cell by cell.
    const amrex::Real scale_factor_mean = 1. / std::sqrt(chi_mean);
    const amrex::Real hubble_mean       = -k_mean / 3.;

    SmallDataIO means_file(m_params.data_path + "means-file", m_dt, m_time,
                           m_restart_time, SmallDataIO::APPEND, m_first_step,
                           ".dat");

    // Drops any data written past this time by a run that was later restarted
    means_file.remove_duplicate_time_data();

    if (m_first_step)
    {
        means_file.write_header_line(
            {"PhiMean", "PiMean", "ScaleFactMean", "HubbleMean", "LapseMean"});
    }

    means_file.write_time_data_line(std::vector<amrex::Real>{
        phi_mean, pi_mean, scale_factor_mean, hubble_mean, lapse_mean});
}

// Bins |field|^2 onto an isotropic k axis and writes the averaged spectrum
inline void
InflationExtraction::print_power_spectrum(const amrex::cMultiFab &field_array,
                                          SmallDataIO &power_spec_file,
                                          const int component)
{
    const int num_bins = params().N / 2;

    // Set up the isotropic k axis bounds
    const amrex::Real kiso_max = std::sqrt(3.) * params().N *
                                 amrex::Math::pi<amrex::Real>() /
                                 params().box_length;
    const amrex::Real dkiso    = std::sqrt(3.) * 2. *
                                 amrex::Math::pi<amrex::Real>() /
                                 params().box_length;

    // Check the stepping along the diagonal is consistent
    if (kiso_max / dkiso - num_bins > InflatonUtils::tolerance)
    {
        amrex::Error("InflationExtraction::print_power_spectrum, "
                     "isotropic k axis is too large.");
    }

    // Set up the isotropic k axis and the binned spectrum
    amrex::Vector<amrex::Real> kiso(num_bins + 1, 0.);
    amrex::Vector<amrex::Real> ps_map(num_bins + 1, 0.);
    amrex::Vector<int> kcount(num_bins + 1, 0);
    for (int s = 0; s <= num_bins; s++)
    {
        kiso[s] = s * dkiso;
    }

    // Device copies of the bins, which the kernel reads and atomically writes
    amrex::Gpu::DeviceVector<amrex::Real> kiso_d(num_bins + 1);
    amrex::Gpu::DeviceVector<amrex::Real> ps_map_d(num_bins + 1);
    amrex::Gpu::DeviceVector<int> kcount_d(num_bins + 1);
    amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, kiso.begin(), kiso.end(),
                          kiso_d.begin());
    amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, ps_map.begin(),
                          ps_map.end(), ps_map_d.begin());
    amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, kcount.begin(),
                          kcount.end(), kcount_d.begin());
    amrex::Gpu::streamSynchronize();

    const amrex::Real *kiso_ptr = kiso_d.data();
    amrex::Real *ps_map_ptr     = ps_map_d.data();
    int *kcount_ptr             = kcount_d.data();

    // Needed to pass the map into the ParallelFor below
    amrex::MFIter::allowMultipleMFIters(true);

    // Local copy so the kernel captures the geometry by value, not via the
    // host `this` pointer
    const InflatonUtils cfg = m_utils;

    const auto &field_arrs = field_array.const_arrays();
    amrex::ParallelFor(
        field_array,
        [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k)
        {
            const amrex::IntVect iv{i, j, k};
            const amrex::Real kmag = cfg.get_kmag(iv);
            const int last_bin     = cfg.m_params.N / 2;

            AMREX_ASSERT_WITH_MESSAGE(
                kmag - kiso_ptr[last_bin] <= InflatonUtils::tolerance,
                "InflationExtraction::print_power_spectrum, |k| above the "
                "isotropic k domain");
            AMREX_ASSERT_WITH_MESSAGE(
                kmag >= kiso_ptr[0],
                "InflationExtraction::print_power_spectrum, |k| below the "
                "isotropic k domain");

            bool binned = false;
            for (int s = 1; s <= last_bin; s++)
            {
                // Bin the power if this mode falls in [kiso[s-1], kiso[s]),
                // or sits exactly on the outermost bin edge
                if ((kmag < kiso_ptr[s] && kmag >= kiso_ptr[s - 1]) ||
                    kmag == kiso_ptr[last_bin])
                {
                    const amrex::Real real_part =
                        field_arrs[box_no](i, j, k, component).real();
                    const amrex::Real imag_part =
                        field_arrs[box_no](i, j, k, component).imag();
                    amrex::Real power =
                        real_part * real_part + imag_part * imag_part;

                    const int bin =
                        (kmag == kiso_ptr[last_bin]) ? last_bin : s - 1;

                    // The r2c transform stores only half of the i range, so
                    // every interior i stands in for its Hermitian partner
                    int count = 1;
                    if (i != 0 && i != last_bin)
                    {
                        power *= 2.;
                        count  = 2;
                    }

                    amrex::Gpu::Atomic::Add(&kcount_ptr[bin], count);
                    amrex::Gpu::Atomic::Add(&ps_map_ptr[bin], power);

                    binned = true;
                    break;
                }
            }

            AMREX_ASSERT_WITH_MESSAGE(
                binned, "InflationExtraction::print_power_spectrum, part of "
                        "the spectrum was not captured by any bin");
            amrex::ignore_unused(binned);
        });

    amrex::Gpu::streamSynchronize();

    // Bring the accumulated bins back to the host
    amrex::Gpu::copyAsync(amrex::Gpu::deviceToHost, ps_map_d.begin(),
                          ps_map_d.end(), ps_map.begin());
    amrex::Gpu::copyAsync(amrex::Gpu::deviceToHost, kcount_d.begin(),
                          kcount_d.end(), kcount.begin());
    amrex::Gpu::streamSynchronize();

    // Each rank has binned only the modes it owns
    amrex::ParallelAllReduce::Sum(kcount.data(),
                                  static_cast<int>(kcount.size()),
                                  amrex::ParallelContext::CommunicatorSub());
    amrex::ParallelAllReduce::Sum(ps_map.data(),
                                  static_cast<int>(ps_map.size()),
                                  amrex::ParallelContext::CommunicatorSub());

    power_spec_file.write_header_line({"power"}, "k");

    for (int s = 0; s < num_bins; s++)
    {
        const amrex::Real avg_power =
            (kcount[s] > 0) ? ps_map[s] / kcount[s] : 0.;
        // |k| is the abscissa of the spectrum, the binned power the datum
        power_spec_file.write_data_line(std::vector<amrex::Real>{avg_power},
                                        (kiso[s] + kiso[s + 1]) / 2.);
    }
}

// Central moment of the requested order, about the supplied mean
inline amrex::Real InflationExtraction::calculate_field_moment_x(
    const amrex::MultiFab &field, const amrex::Real mean, const int moment,
    const int component) const
{
    const amrex::Real vol = std::pow(params().N, 3.);

    amrex::ReduceOps<amrex::ReduceOpSum> reduce_op;
    amrex::ReduceData<amrex::Real> reduce_data(reduce_op);
    using ReduceTuple = typename decltype(reduce_data)::Type;

    for (amrex::MFIter mfi(field); mfi.isValid(); ++mfi)
    {
        const amrex::Box &box = mfi.validbox();
        auto const &field_arr = field.const_array(mfi);
        reduce_op.eval(
            box, reduce_data,
            [=] AMREX_GPU_DEVICE(int i, int j, int k) -> ReduceTuple
            {
                return {std::pow(field_arr(i, j, k, component) - mean, moment)};
            });
    }

    // reduce_data.value() is local to this rank
    amrex::Real sum = amrex::get<0>(reduce_data.value());
    amrex::ParallelAllReduce::Sum(sum,
                                  amrex::ParallelContext::CommunicatorSub());

    if (sum == 0.)
    {
        return 0.;
    }
    // The second moment is reported as a standard deviation
    return (moment == 2) ? std::sqrt(sum / vol) : sum / vol;
}

inline amrex::Vector<amrex::Real> InflationExtraction::print_moment(
    const amrex::MultiFab &field, const amrex::Vector<std::string> &names,
    const amrex::Vector<int> &moment_orders, SmallDataIO &file)
{
    const int num_comps = field.nComp();
    amrex::Vector<amrex::Real> stdev(num_comps, 0.);

    if (moment_orders.empty())
    {
        return stdev;
    }

    // check_params() guarantees the orders are ascending and within 1..4
    const int max_order   = moment_orders.back();
    const amrex::Real vol = std::pow(params().N, 3.);

    const auto start     = moment_orders.begin();
    const auto mean_itr  = std::find(start, moment_orders.end(), 1);
    const auto stdev_itr = std::find(start, moment_orders.end(), 2);
    const auto skew_itr  = std::find(start, moment_orders.end(), 3);
    const auto kurt_itr  = std::find(start, moment_orders.end(), 4);

    const int num_entries = num_comps * static_cast<int>(moment_orders.size());
    amrex::Vector<amrex::Real> data_to_print(num_entries, 0.);
    amrex::Vector<std::string> headers(num_entries, "");

    for (int comp = 0; comp < num_comps; comp++)
    {
        // MultiFab::sum is already collective across ranks
        const amrex::Real mean = field.sum(comp) / vol;
        if (mean_itr != moment_orders.end())
        {
            assign_statistics_data(headers, names[comp] + "-mean",
                                   data_to_print, mean, comp, num_comps,
                                   mean_itr, start, m_first_step);
        }

        if (max_order >= 2)
        {
            stdev[comp] = calculate_field_moment_x(field, mean, 2, comp);
            if (stdev_itr != moment_orders.end())
            {
                assign_statistics_data(
                    headers, names[comp] + "-stdev", data_to_print, stdev[comp],
                    comp, num_comps, stdev_itr, start, m_first_step);
            }
        }

        // Skewness and kurtosis are normalised by the standard deviation, so
        // they are only meaningful for a field with non-zero variance
        const bool has_spread = (stdev[comp] != 0.);

        if (max_order >= 3 && skew_itr != moment_orders.end())
        {
            const amrex::Real skew =
                has_spread ? calculate_field_moment_x(field, mean, 3, comp) /
                                 std::pow(stdev[comp], 3.)
                           : 0.;
            assign_statistics_data(headers, names[comp] + "-skew",
                                   data_to_print, skew, comp, num_comps,
                                   skew_itr, start, m_first_step);
        }

        if (max_order >= 4 && kurt_itr != moment_orders.end())
        {
            const amrex::Real kurt =
                has_spread ? calculate_field_moment_x(field, mean, 4, comp) /
                                 std::pow(stdev[comp], 4.)
                           : 0.;
            assign_statistics_data(headers, names[comp] + "-kurt",
                                   data_to_print, kurt, comp, num_comps,
                                   kurt_itr, start, m_first_step);
        }
    }

    if (m_first_step)
    {
        file.write_header_line(headers);
    }
    file.write_time_data_line(data_to_print);

    return stdev;
}

/* Main routine */

inline void InflationExtraction::extract(const amrex::MultiFab &state)
{
    BL_PROFILE("InflationExtraction::extract");

    FilesystemTools::ensure_directory_exists(m_params.data_path);

    // The background averages are read straight off the state, so they are
    // taken before the curvature and polarisation fields are reconstructed
    if (m_params.calc_background_means != 0)
    {
        extract_background_means(state);
    }

    const int step = static_cast<int>(std::round(m_time / m_dt));

    const bool want_spectra    = (m_params.calc_binned_power_spectrum != 0) &&
                                 (step % m_params.spec_interval == 0);
    const bool want_statistics = (m_params.calc_higher_order_statistics != 0) &&
                                 !m_params.moment_orders.empty();

    if (!want_spectra && !want_statistics)
    {
        return;
    }

    // Configuration-space fields, and the Fourier-space fields behind them
    // when a spectrum is wanted this step
    const amrex::BoxArray &sba            = state.boxArray();
    const amrex::DistributionMapping &sdm = state.DistributionMap();
    amrex::MultiFab hs_x(sba, sdm, 2, 0);
    amrex::MultiFab R_x(sba, sdm, 1, 0);
    hs_x.setVal(0.0);
    R_x.setVal(0.0);

    amrex::cMultiFab hs_k;
    amrex::cMultiFab R_k;

    DerivedVariables extractor;
    extractor.extract_hs_and_R(hs_x, R_x, state, want_spectra ? &hs_k : nullptr,
                               want_spectra ? &R_k : nullptr);

    if (want_spectra)
    {
        const std::string spec_path = make_subdirectory("spectra");

        for (int comp = 0; comp < hs_k.nComp(); comp++)
        {
            // hs_k holds (hplus, hcross); var_names leads with R
            SmallDataIO spectrum_file(
                spec_path + "spectrum-" + var_names[comp + 1] + "-time-", m_dt,
                m_time, m_restart_time, SmallDataIO::NEW, m_first_step, ".dat");
            print_power_spectrum(hs_k, spectrum_file, comp);
        }

        SmallDataIO spectrum_file(
            spec_path + "spectrum-" + var_names[0] + "-time-", m_dt, m_time,
            m_restart_time, SmallDataIO::NEW, m_first_step, ".dat");
        print_power_spectrum(R_k, spectrum_file, 0);
    }

    if (!want_statistics)
    {
        return;
    }

    // Gather R and the polarisations into one MultiFab, in var_names order
    const int output_comps = R_x.nComp() + hs_x.nComp();
    amrex::MultiFab out_mf(sba, sdm, output_comps, 0);
    amrex::MultiFab::Copy(out_mf, R_x, 0, 0, R_x.nComp(), 0);
    amrex::MultiFab::Copy(out_mf, hs_x, 0, R_x.nComp(), hs_x.nComp(), 0);

    SmallDataIO stats_file(m_params.data_path + "field-statistics", m_dt,
                           m_time, m_restart_time, SmallDataIO::APPEND,
                           m_first_step, ".dat");
    stats_file.remove_duplicate_time_data();

    const amrex::Vector<amrex::Real> stdevs =
        print_moment(out_mf, var_names, m_params.moment_orders, stats_file);

    // The tensor-to-scalar ratio is built from the standard deviations, so it
    // is only available when the second moment was requested
    const bool have_stdev =
        std::find(m_params.moment_orders.begin(), m_params.moment_orders.end(),
                  2) != m_params.moment_orders.end();

    if (have_stdev && stdevs[0] != 0.)
    {
        SmallDataIO ts_file(m_params.data_path + "tensor-scalar-ratio", m_dt,
                            m_time, m_restart_time, SmallDataIO::APPEND,
                            m_first_step, ".dat");
        ts_file.remove_duplicate_time_data();

        if (m_first_step)
        {
            ts_file.write_header_line(
                {"T/S ratio (plus)", "T/S ratio (cross)"});
        }
        ts_file.write_time_data_line(std::vector<amrex::Real>{
            stdevs[1] / stdevs[0], stdevs[2] / stdevs[0]});
    }
}

#endif /* INFLATIONEXTRACTION_IMPL_HPP_ */

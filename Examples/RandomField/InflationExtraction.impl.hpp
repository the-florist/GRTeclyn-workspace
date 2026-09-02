/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef INFLATIONEXTRACTION_IMPL_HPP_
#define INFLATIONEXTRACTION_IMPL_HPP_

void InflationExtraction::set_up(int a_state_index)
{
    int num_ghosts = 0;

    auto &derive_lst     = amrex::AmrLevel::get_derive_lst();
    const auto &desc_lst = amrex::AmrLevel::get_desc_lst();

    const auto &field_names = InflationExtraction::var_names;

    // Add Constraints to the derive list
    derive_lst.add(
        "InflationFields", amrex::IndexType::TheCellType(),
        static_cast<int>(field_names.size()), field_names,
        InflationExtraction::compute_mf, [=](const amrex::Box &box)
        { return amrex::grow(box, num_ghosts); }, &amrex::cell_quartic_interp);

    derive_lst.addComponent("InflationFields", desc_lst, a_state_index, 0,
                            NUM_VARS);
}

/* Helper functions */

// Makes subdirectories in data/
inline std::string 
InflationExtraction::make_subdirectory(const std::string base, 
                                       const std::string dir, 
                                       const int is_first_step) const
{
    std::string new_path = base+"../"+dir+"/";
    if(is_first_step)
    {
        if (!FilesystemTools::directory_exists(base))
        {
            amrex::Print() << "InflationExtraction::make_subdirectory, "
                              "Directory creation failed for "; 
            amrex::Print() << new_path << "\n";
        }
        else if (!FilesystemTools::directory_exists(new_path))
        {
            FilesystemTools::mkdir_recursive(new_path);
        }
    }
    return new_path;
}

// Creates a custom data file layout 
inline void
InflationExtraction::assign_statistics_data(
    amrex::Vector<std::string> &header_storage, const std::string name,
    amrex::Vector<amrex::Real> &data_storage, const amrex::Vector<amrex::Real> data,
    const int component, const int num_comps,
    const amrex::Vector<int>::const_iterator itr,
    const amrex::Vector<int>::const_iterator start, const int is_first_step)
{
    int loc = component + num_comps*(itr - start);
    if(is_first_step)
    {
        header_storage[loc] =  name; 
    }
    data_storage[loc] = data[component];
}

// Calculates and prints the power spectrum
inline void InflationExtraction::print_power_spectrum(
    const amrex::cMultiFab &field_array, SmallDataIO &power_spec_file,
    const int component)
{ 
    // Set up the isotropic k axis bounds
    amrex::Real kiso_max = std::sqrt(3.) * inflt_methods.N *
                           amrex::Math::pi<amrex::Real>() / inflt_methods.L;
    amrex::Real dkiso = sqrt(3.) * 2. * amrex::Math::pi<amrex::Real>() / inflt_methods.L;

    // check the stepping along the diagonal is consistent
    if (kiso_max/dkiso - inflt_methods.N/2 > InflationUtils::tolerance)
    {
        amrex::Error("InflationExtraction::print_power_spectrum "
                     "Isotropic k axis is too large.");
    }

    // Set up isotropic k axis and PS map
    amrex::Vector<amrex::Real> kiso(inflt_methods.N / 2 + 1, 0.);
    amrex::Vector<amrex::Real> ps_map(inflt_methods.N/2 + 1, 0.);
    amrex::Vector<int> kcount(inflt_methods.N/2 + 1, 0);
    for (int s=0; s<=inflt_methods.N/2; s++) { kiso[s] = s*dkiso; }

    // Device copies of the bins, which the kernel reads and atomically writes
    amrex::Gpu::DeviceVector<amrex::Real> kiso_d(inflt_methods.N/2 + 1);
    amrex::Gpu::DeviceVector<amrex::Real> ps_map_d(inflt_methods.N/2 + 1);
    amrex::Gpu::DeviceVector<int> kcount_d(inflt_methods.N/2 + 1);
    amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, kiso.begin(), kiso.end(),
                          kiso_d.begin());
    amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, ps_map.begin(), ps_map.end(),
                          ps_map_d.begin());
    amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, kcount.begin(), kcount.end(),
                          kcount_d.begin());
    amrex::Gpu::streamSynchronize();
    const amrex::Real* kiso_ptr = kiso_d.data();
    amrex::Real* ps_map_ptr = ps_map_d.data();
    int* kcount_ptr = kcount_d.data();

    // Needed to pass the map to the amrex::ParallelFor loop
    amrex::MFIter::allowMultipleMFIters(true);

    // Slice to the POD base so the kernel captures config by value
    const InflationConfig cfg = inflt_methods;

    // Loop to bin the power spectrum at each point
    const auto& field_arrs = field_array.arrays();
    amrex::ParallelFor(field_array, [=]
                AMREX_GPU_DEVICE (int bx, int i, int j, int k)
        {
            // Check to see if you're in a ghost cell
            amrex::IntVect iv{i, j, k};
            amrex::Real kmag = cfg.get_kmag(iv);

            // make sure you're still in the domain
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
                kmag - kiso_max <= InflationUtils::tolerance,
                "InflationExtraction::print_power_spectrum, magnitude "
                "larger than (N/2,N/2,N/2)");

            // Loop over the isotropic axis
            for (int s=1; s<=cfg.N/2; s++)
            {
                // If smaller than the smallest bin
                AMREX_ALWAYS_ASSERT_WITH_MESSAGE(kmag >= kiso_ptr[0],
                        "InflationExtraction::print_power_spectrum, "
                        "kmag below the kiso domain");

                // If you're larger than the largest bin
                AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
                    kmag - kiso_ptr[cfg.N/2] <= InflationUtils::tolerance,
                        "InflationExtraction::print_power_spectrum, "
                        "kmag above the kiso domain");

                // If you're somewhere in the middle, bin the power here
                if ((kmag < kiso_ptr[s] && kmag >= kiso_ptr[(s-1)])
                        || kmag == kiso_ptr[cfg.N/2])
                {
                    amrex::Real power =
                        (std::pow(field_arrs[bx](i, j, k, component).real(), 2.0)
                         + std::pow(field_arrs[bx](i, j, k, component).imag(), 2.0));

                    int comp = (kmag == kiso_ptr[cfg.N/2]) ? cfg.N/2 : s - 1;

                    int count = 1;
                    if (i != 0 && i != cfg.N/2) { power *= 2.; count = 2; }

                    amrex::Gpu::Atomic::Add(&kcount_ptr[comp], count);
                    amrex::Gpu::Atomic::Add(&ps_map_ptr[comp], power);

                    break;
                }

                // If you've reached the largest bin but not been captured
                AMREX_ALWAYS_ASSERT_WITH_MESSAGE(s <= cfg.N/2,
                        "InflationExtraction::print_power_spectrum, "
                        "part of the spectrum isn't captured");
            }
        });

    amrex::Gpu::streamSynchronize();

    // Bring the accumulated bins back to the host
    amrex::Gpu::copyAsync(amrex::Gpu::deviceToHost, ps_map_d.begin(),
                          ps_map_d.end(), ps_map.begin());
    amrex::Gpu::copyAsync(amrex::Gpu::deviceToHost, kcount_d.begin(),
                          kcount_d.end(), kcount.begin());
    amrex::Gpu::streamSynchronize();

    amrex::ParallelAllReduce::Sum(kcount.data(), static_cast<int>(kcount.size()),
                                  amrex::ParallelContext::CommunicatorSub());
    amrex::ParallelAllReduce::Sum(ps_map.data(), static_cast<int>(ps_map.size()), 
                                  amrex::ParallelContext::CommunicatorSub());

    for(int s = 0; s < inflt_methods.N/2; s++)
    {
        const amrex::Real avg_power = (kcount[s] > 0) ? ps_map[s] / kcount[s] : 0.;
        const std::vector<amrex::Real> data{(kiso[s] + kiso[s+1]) / 2., avg_power};
        power_spec_file.write_data_line(data);
    }
}

// Finds statistical moment x of given MultiFab
inline amrex::Real 
InflationExtraction::calculate_field_moment_x(const amrex::MultiFab &field, 
                                              const amrex::Vector<amrex::Real> mean, 
                                              const int moment, const int component)
{
    const amrex::Real vol = std::pow(inflt_methods.N, 3.);
    const amrex::Real mean_comp = mean[component];

    amrex::ReduceOps<amrex::ReduceOpSum> reduce_op;
    amrex::ReduceData<amrex::Real> reduce_data(reduce_op);
    using ReduceTuple = typename decltype(reduce_data)::Type;

    for (amrex::MFIter mfi(field); mfi.isValid(); ++mfi)
    {
        const amrex::Box &bx = mfi.validbox();
        auto const &arr = field.const_array(mfi);
        reduce_op.eval(bx, reduce_data,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
            {
                return { std::pow(arr(i, j, k, component) - mean_comp, moment) };
            });
    }

    amrex::Real sum = amrex::get<0>(reduce_data.value());
    amrex::ParallelAllReduce::Sum(sum, amrex::ParallelContext::CommunicatorSub());

    // Normalise and return moment x
    if (sum == 0) { return 0; }
    else if(moment == 2) { return sqrt(sum/vol); }
    else { return sum/vol; }
}

// Calculates and prints requested moments (any between 1 and 4)
inline amrex::Vector<amrex::Real> 
InflationExtraction::print_moment(amrex::MultiFab &field, 
                                  const amrex::Vector<std::string> names,  
                                  const amrex::Vector<int> &moment_orders, 
                                  SmallDataIO &file, 
                                  const int is_first_step)
{
    // Trap instance where the user requests too large a moment
    for(const auto moment : moment_orders)
    {
        if(moment > 4) 
        { 
            amrex::Error("InflationExtraction::print_moment "
                  "Chosen moment order has not been implemented");
        }
    }

    // Allocate arrays to store each moment
    const int nc = field.nComp();
    const amrex::Real vol = std::pow(inflt_methods.N, 3.);
    amrex::Vector<amrex::Real> means(nc, 0.);
    amrex::Vector<amrex::Real> stdev(nc, 0.);
    amrex::Vector<amrex::Real> skew(nc, 0.);
    amrex::Vector<amrex::Real> kurt(nc, 0.);

    // Find iterators, which determine which moments are requested and their ordering
    amrex::Vector<int>::const_iterator start = moment_orders.begin();
    amrex::Vector<int>::const_iterator mean_itr = std::find(moment_orders.begin(), 
                                                            moment_orders.end(), 1);
    amrex::Vector<int>::const_iterator stdev_itr = std::find(moment_orders.begin(), 
                                                             moment_orders.end(), 2);
    amrex::Vector<int>::const_iterator skew_itr = std::find(moment_orders.begin(), 
                                                            moment_orders.end(), 3);
    amrex::Vector<int>::const_iterator kurt_itr = std::find(moment_orders.begin(), 
                                                            moment_orders.end(), 4);

    // Allocate vectors to store header line and data lines
    amrex::Vector<amrex::Real> data_to_print(nc * moment_orders.size(), 0.);
    amrex::Vector<std::string> headers(nc * moment_orders.size(), "");

    for (int comp = 0; comp < nc; comp++)
    {
        means[comp] = field.sum(comp)/vol;
        if(mean_itr != moment_orders.end())
        {
            assign_statistics_data(headers, names[comp]+"-mean", data_to_print,
                                    means, comp, nc, mean_itr, start, is_first_step);
        }

        if(moment_orders.back() != 1)
        {
            stdev[comp] = calculate_field_moment_x(field, means, 2, comp);
            if(stdev_itr != moment_orders.end())
            {
                assign_statistics_data(headers, names[comp]+"-stdev", data_to_print,
                                        stdev, comp, nc, stdev_itr, start, is_first_step);
            }

            if(moment_orders.back() != 2)
            {
                skew[comp] = calculate_field_moment_x(field, means, 3, comp);
                skew[comp] /= std::pow(stdev[comp], 3.);

                if(skew_itr != moment_orders.end())
                {
                    assign_statistics_data(headers, names[comp]+"-skew",
                                           data_to_print, skew, comp, nc,
                                           skew_itr, start, is_first_step);
                }

                if(moment_orders.back() != 3)
                {
                    kurt[comp] = calculate_field_moment_x(field, means, 4, comp);
                    kurt[comp] /= std::pow(stdev[comp], 4.);

                    assign_statistics_data(headers, names[comp]+"-kurt",
                                           data_to_print, kurt, comp, nc,
                                           kurt_itr, start, is_first_step);
                }
            }
        }
    }

    // Write header and data lines
    if(is_first_step) 
    { 
        file.write_header_line(headers); 
    }

    file.write_time_data_line(data_to_print);
    
    return stdev;
}

/* Main functions */

// Extract R and hs in configuration space from the BSSN variables
inline void InflationExtraction::extract_hs_and_R(amrex::MultiFab &hs, 
                                                  amrex::MultiFab &R, 
                                                  const amrex::MultiFab &state, 
                                                  const bool print_spec = false)
{
    // Extract amrex::MultiFab ingredients from state
    amrex::BoxArray sba = state.boxArray();
    amrex::DistributionMapping sdm = state.DistributionMap();
    if (sba != hs.boxArray() 
        || sdm != hs.DistributionMap())
    {
        amrex::Error("InflationExtraction::extract_hs_and_R "
              "source and output BA or SDM do not match");
    }

    // 0: scalar field
    // 1: conformal factor
    amrex::MultiFab scalars_x(sba, sdm, 2, 0);
    amrex::MultiFab gij_x(sba, sdm, 6, 0);

    // Copy the spatial metric from the state
    Copy(gij_x, state, c_h11, InflationUtils::lut[0][0], 1, 0);
    Copy(gij_x, state, c_h12, InflationUtils::lut[0][1], 1, 0);
    Copy(gij_x, state, c_h13, InflationUtils::lut[0][2], 1, 0);
    Copy(gij_x, state, c_h22, InflationUtils::lut[1][1], 1, 0);
    Copy(gij_x, state, c_h23, InflationUtils::lut[1][2], 1, 0);
    Copy(gij_x, state, c_h33, InflationUtils::lut[2][2], 1, 0);

    constexpr int m_c_phi = 0;
    constexpr int m_c_chi = 1;
    Copy(scalars_x, state, c_phi, m_c_phi, 1, 0);
    Copy(scalars_x, state, c_chi, m_c_chi, 1, 0);

    // Find background quantities needed to extract \cal R
    const int vol = std::pow(inflt_methods.N, 3);
    const amrex::Real K_bar = state.sum(c_K)/vol;
    const amrex::Real alpha_bar = state.sum(c_lapse)/vol;
    const amrex::Real Pi_bar = state.sum(c_Pi)/vol;
    const amrex::Real phi_bar = state.sum(c_phi)/vol;
    const amrex::Real chi_bar = state.sum(c_chi)/vol;

    // Remove background from scalar field
    scalars_x.plus(-phi_bar, m_c_phi, 1);
    scalars_x.plus(-chi_bar, m_c_chi, 1);
    scalars_x.mult(1./inflt_methods.norm());

    // Undo the normalisation and BSSN-CPT conversion
    for (int l=0; l<3; l++) { gij_x.plus(-1., InflationUtils::lut[l][l], 1); }
    gij_x.mult(1./inflt_methods.norm());

    // Set up the problem domain in Fourier space
    // And impose that MPI ranks only slice along the i index (for Nyquist conditions)
    amrex::IntVect domain_low(0, 0, 0);
    amrex::IntVect k_domain_high(inflt_methods.N/2, inflt_methods.N-1, inflt_methods.N-1);
    amrex::Box k_domain(domain_low, k_domain_high);
    constexpr amrex::Array<bool, AMREX_SPACEDIM> slicing{true, false, false};
    amrex::BoxArray kba =
        decompose(k_domain, amrex::ParallelContext::NProcsAll(), slicing);
    amrex::DistributionMapping kdm(kba);

    // Set up the arrays to store the Fourier data sets
    amrex::cMultiFab hs_k(kba, kdm, 2, 0);
    amrex::cMultiFab gij_k(kba, kdm, 6, 0);
    amrex::cMultiFab scalars_k(kba, kdm, 2, 0);
    amrex::cMultiFab R_k(kba, kdm, 1, 0);

    hs_k.setVal(0.0);
    gij_k.setVal(0.0);
    scalars_k.setVal(0.0);
    R_k.setVal(0.0);

    // Set up the FFT
    amrex::IntVect x_domain_high(inflt_methods.N-1, inflt_methods.N-1, inflt_methods.N-1);
    amrex::Box x_domain(domain_low, x_domain_high);
    amrex::FFT::R2C<amrex::Real> tensor_fft(
        x_domain, amrex::FFT::Info().setBatchSize(gij_k.nComp()));
    amrex::FFT::R2C<amrex::Real> scalar_fft(
        x_domain, amrex::FFT::Info().setBatchSize(scalars_k.nComp()));
    amrex::FFT::R2C<amrex::Real> mode_fn_fft(
        x_domain, amrex::FFT::Info().setBatchSize(hs_k.nComp()));
    amrex::FFT::R2C<amrex::Real> R_fft(
        x_domain, amrex::FFT::Info().setBatchSize(R_k.nComp()));

    // Perform the fft
    tensor_fft.forward(gij_x, gij_k);
    scalar_fft.forward(scalars_x, scalars_k);

    // Normalise the fft (fftw style)
    for(int comp = 0; comp < 6; comp++)
    {
        gij_k.mult(std::pow(inflt_methods.N, -3.), comp, 1);
    }
    for(int comp = 0; comp < 2; comp++)
    {
        scalars_k.mult(std::pow(inflt_methods.N, -3.), comp, 1);
    }

    const auto& hs_arrs = hs_k.arrays();
    const auto& gij_arrs = gij_k.arrays();
    const auto& scalars_arrs = scalars_k.arrays();
    const auto& R_k_arrs = R_k.arrays();

    // Slice to the POD base so the kernel captures config by value
    const InflationConfig cfg = inflt_methods;

    amrex::ParallelFor(gij_k, [=]
                AMREX_GPU_DEVICE (int bx, int i, int j, int k)
        {
            amrex::IntVect iv{i, j, k};
            amrex::Real kmag = cfg.get_kmag(iv);

            if (iv != amrex::IntVect{0, 0, 0})
            {
                const auto [eplus, ecross] = cfg.calculate_polarisation_tensors(iv);

                // Find basis tensors and do the Fourier trick
                for (int l=0; l<3; l++) for (int p=0; p<3; p++)
                {
                    hs_arrs[bx](i, j, k, 0) +=
                        (gij_arrs[bx](i, j, k, InflationUtils::lut[l][p])
                         * eplus(l, p))/2.;
                    hs_arrs[bx](i, j, k, 1) +=
                        (gij_arrs[bx](i, j, k, InflationUtils::lut[l][p])
                         * ecross(l, p))/2.;
                }

                // Calculate the TT and scalar-(vector) components of the 
                // metric, by reconstructing hij and subtracting it from \tilde{gamma}_ij
                Tensor<2, amrex::GpuComplex<amrex::Real>> hij, hSV;
                for (int l=0; l<3; l++) for (int p=0; p<3; p++)
                {
                    hij[l][p] = (hs_arrs[bx](i, j, k, 0) * eplus(l, p)
                                + hs_arrs[bx](i, j, k, 1) * ecross(l, p));
                    hSV[l][p] =
                        gij_arrs[bx](i, j, k, InflationUtils::lut[l][p]) - hij[l][p];
                }

                // Extract R according to the scheme detailed in
                // Appendix B (Eq. B1) of arxiv:2502.06783, using hSV as the 
                // spatial metric instead of \tilde{gamma}_ij
                if(cfg.scalar_init)
                {
                    // Find the unitful k vector
                    amrex::GpuArray<amrex::Real, 3> iv_k{
                        static_cast<amrex::Real>(iv[0]),
                        static_cast<amrex::Real>(iv[1]),
                        static_cast<amrex::Real>(iv[2])};
                    iv_k[1] = cfg.invert_index_with_sign(iv_k[1]);
                    iv_k[2] = cfg.invert_index_with_sign(iv_k[2]);

                    for(auto& k_comp : iv_k)
                    {
                        k_comp *= 2. * amrex::Math::pi<amrex::Real>() / cfg.L;
                    }
                    amrex::GpuComplex<amrex::Real> Phi = 0;

                    // converstion from chi and gamma_ij -> Phi
                    for(int l=0; l<3; l++) for(int p=0; p<3; p++)
                    {
                        Phi += (iv_k[l] * iv_k[p] * hSV[l][p])/std::pow(kmag, 2.);
                    }
                    Phi *= 1./4.;
                    Phi += 0.5 * (scalars_arrs[bx](i, j, k, m_c_chi));

                    // Combine the above to find R(k)
                    R_k_arrs[bx](i, j, k, 0) =
                        Phi - ((K_bar/3.) * scalars_arrs[bx](i, j, k, m_c_phi)
                               / alpha_bar / Pi_bar);
                }
            }
        });

    amrex::Gpu::streamSynchronize();

    // Prepare to IFT the polarisation fields and R field
    inflt_methods.apply_nyquist_conditions(hs_k);
    inflt_methods.apply_nyquist_conditions(R_k);

    // Find the binned PS for each mode function and print to data/
    GRParmParse extraction_pp("extraction");
    int spec_interval = 100;
    extraction_pp.query("spec_interval", spec_interval);
    if ((print_spec) && (static_cast<int>(time/dt)
                         % spec_interval == 0))
    {
        amrex::Print() << "Time step at print: "
                       << static_cast<int>(std::round(time/dt)) << "\n";
        std::string spec_path = make_subdirectory(m_data_path, "spectra", first_step);

        for(int comp = 0; comp < hs_k.nComp(); comp++)
        {
            SmallDataIO spectrum_file(
                spec_path+"spectrum-comp-"+std::to_string(comp)+"-time-",
                dt, time, restart_time, SmallDataIO::NEW, first_step, ".dat");
            print_power_spectrum(hs_k, spectrum_file, comp);
        }

        SmallDataIO spectrum_file(spec_path+"spectrum-Rk-time-", dt, time, 
                                  restart_time, SmallDataIO::NEW, first_step, ".dat");
        print_power_spectrum(R_k, spectrum_file, 0);
    }

    // Fourier transform
    mode_fn_fft.backward(hs_k, hs);
    R_fft.backward(R_k, R);

    // Confirm Parseval's theorem holds between config and Fourier space
    // (before applying the physical normalisation)
    if(print_spec)
    {
        TensorTests::Test_Parsevals_thm(hs, hs_k, inflt_methods.N);
        TensorTests::Test_Parsevals_thm(R, R_k, inflt_methods.N);
    }

    // Apply physical normalisation
    hs.mult(inflt_methods.norm());
    R.mult(inflt_methods.norm());
}

// Put R and hs into plotfiles
// DeriveFuncMF callback: build an extractor and fill the plotfile output with
// R, hplus, hcross. src_mf arrives already FillPatch-ed by the framework.
inline void InflationExtraction::compute_mf(
    amrex::MultiFab &out_mf, int dcomp, int /*ncomp*/,
    const amrex::MultiFab &src_mf, const amrex::Geometry & /*geomdata*/,
    amrex::Real /*time*/, const int * /*bcrec*/, int /*level*/)
{
    InflationExtraction extractor;
    extractor.derive(src_mf, out_mf, dcomp);
}

inline void InflationExtraction::derive(const amrex::MultiFab &source,
                                        amrex::MultiFab &out, const int dcomp)
{
    BL_PROFILE("InflationExtraction::derive");

    // Make a multifab to store config space mode functions
    amrex::BoxArray oba = out.boxArray();
    amrex::DistributionMapping odm = out.DistributionMap();
    amrex::MultiFab hs_x(oba, odm, 2, 0);
    amrex::MultiFab R_x(oba, odm, 1, 0);
    hs_x.setVal(0.0);
    R_x.setVal(0.0);

    // print_spec = false: no data-file side effects on the plotfile path
    extract_hs_and_R(hs_x, R_x, source, false);

    const auto& hs_arrs = hs_x.arrays();
    const auto& R_arrs = R_x.arrays();
    const auto& out_arrs = out.arrays();

    amrex::ParallelFor(hs_x, [=, this] AMREX_GPU_DEVICE
                (int bx, int i, int j, int k)
        {
            const amrex::IntVect iv{i, j, k};
            out_arrs[bx](iv, dcomp) = R_arrs[bx](i, j, k);
            out_arrs[bx](iv, dcomp + 1) = hs_arrs[bx](i, j, k, 0);
            out_arrs[bx](iv, dcomp + 2) = hs_arrs[bx](i, j, k, 1);
        });

    amrex::Gpu::streamSynchronize();
}

// Find spectrum and higher-order statistics on R and hs
inline void InflationExtraction::extract(const amrex::MultiFab &state)
{
    BL_PROFILE("InflationExtraction::extract");

    // Make a multifab to store config space mode functions
    amrex::BoxArray sba = state.boxArray();
    amrex::DistributionMapping sdm = state.DistributionMap();
    amrex::MultiFab hs_x(sba, sdm, 2, 0);
    amrex::MultiFab R_x(sba, sdm, 1, 0);
    hs_x.setVal(0.0);
    R_x.setVal(0.0);

    extract_hs_and_R(hs_x, R_x, state, true);

    // Find mode functions in configuration space if requested,
    // and find the statistics (orders 1-4) of the polarisation fields 
    // and the R field. 
    // And print the tensor-to-scalar ratio if requested.
    const int output_comps = hs_x.nComp() + R_x.nComp();
    amrex::MultiFab out_MF(hs_x.boxArray(), hs_x.DistributionMap(), output_comps, 0);
    amrex::MultiFab::Copy(out_MF, R_x, 0, 0, R_x.nComp(), 0);
    amrex::MultiFab::Copy(out_MF, hs_x, 0, R_x.nComp(), hs_x.nComp(), 0);

    // Read the requested moment orders
    GRParmParse extraction_pp("extraction");
    amrex::Vector<int> orders;
    extraction_pp.queryarr("moments_to_print", orders);

    // Calculate and print field moments
    amrex::Vector<amrex::Real> stdevs;
    SmallDataIO stats_file(m_data_path+"field-statistics", dt, time,
                            restart_time, SmallDataIO::APPEND, first_step, ".dat");

    if (!orders.empty())
    {
        stdevs = print_moment(out_MF, var_names, orders,
                              stats_file, first_step);
    }

    // Calculate and print tensor to scalar ratio (integrated PS)
    if (std::find(orders.begin(), orders.end(), 2)
        != orders.end())
    {
        SmallDataIO ts_file(m_data_path+"tensor-scalar-ratio", dt, time,
                            restart_time, SmallDataIO::APPEND, first_step, ".dat");

        if(first_step)
        {
            ts_file.write_header_line({"T/S ratio (plus)", "T/S ratio (cross)"});
        }

        ts_file.write_time_data_line({stdevs[1] / stdevs[0], stdevs[2] / stdevs[0]});
    }
}

#endif /* INFLATIONEXTRACTION_IMPL_HPP_ */
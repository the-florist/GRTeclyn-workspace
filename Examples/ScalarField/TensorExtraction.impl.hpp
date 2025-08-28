/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#if !defined(TENSOREXTRACTION_HPP_)
#error "This file should only be included via TensorExtraction.hpp"
#endif

#ifndef TENSOREXTRACTION_IMPL_HPP_
#define TENSOREXTRACTION_IMPL_HPP_

// Calculates and prints the power spectrum
template <class matter_t>
inline void TensorExtraction<matter_t>::print_power_spectrum(cMultiFab &field_array, SmallDataIO &power_spec_file, const int component)
{ 
    // Set up the isotropic k axis bounds
    double kiso_max = std::sqrt(3.) * N * M_PI / L;
    double dkiso = sqrt(3.)*2.*M_PI/L;
    double tolerance = 1.e-12;

    // check the stepping along the diagonal is consistent
    if (kiso_max/dkiso - N/2 > tolerance)
    {
        Error("TensorExtraction::print_power_spectrum Isotropic k axis is too large.");
    }
    // check you aren't sampling above the max sampling rate
    else if (m_params.bin_number > kiso_max/dkiso)
    {
        Error("TensorExtraction::print_power_spectrum Bin number is too large.");
    }

    // Set up isotropic k axis and PS map
    double dk_to_bin = (double)m_params.bin_number/((double)N/2);
    double kmag = 0.;
    Vector<Real> kiso(N/2+1, 0.);

    Vector<Real> ps_map(m_params.bin_number+1, 0.);
    Vector<int> kcount(m_params.bin_number+1, 0);

    for (int s=0; s<=N/2; s++) { kiso[s] = s*dkiso; }

    // Loop to bin the power spectrum at each point
    MFIter::allowMultipleMFIters(true); // Needed to pass the map to the ParallelFor loop
    for (MFIter mfi(field_array); mfi.isValid(); ++mfi) 
    {
        // Make a pointer to the mode function at this MF box
        Array4<GpuComplex<Real>> const& field_ptr = field_array.array(mfi);
        const Box& bx = mfi.fabbox();

        amrex::ParallelFor(bx, [=, &ps_map, &kcount] AMREX_GPU_DEVICE (int i, int J, int K) noexcept
        {
            IntVect iv{i, J, K};

            // Check to see if you're in a ghost cell
            bool in_ghost_index = is_ghost_index(iv);
            if(!in_ghost_index)
            {
                int j = invert_index(J);
                int k = invert_index(K);
                double kmag = std::sqrt(i*i + j*j + k*k) * 2 * M_PI / L;

                // make sure you're still in the domain
                if(kmag > kiso_max) 
                { 
                    Print() << iv << "\n";
                    Print() << kmag << "," << kiso_max << "\n";
                    Error("TensorExtraction::print_power_spectrum Found magnitude larger than (N/2,N/2,N/2)."); 
                }

                // Loop over the isotropic axis
                for (int s=1; s<=N/2; s++) 
                {
                    // If smaller than the smallest bin
                    if(kmag < kiso[0])
                    {
                        Print() << iv << "\n";
                        Error("TensorExtraction::print_power_spectrum kmag below the kiso domain.");
                    }

                    // If you're larger than the largest bin
                    else if(kmag - kiso[N/2] > tolerance)
                    {
                        Print() << iv << "\n";
                        Error("TensorExtraction::print_power_spectrum kmag above the kiso domain.");
                    }

                    // If you're somewhere in the middle
                    else if (kmag < kiso[s] && kmag >= kiso[(s-1)]) 
                    {
                        Real power = (std::pow(field_ptr(i, j, k, component).real(), 2.0) 
                                    + std::pow(field_ptr(i, j, k, component).imag(), 2.0));
                        
                        Gpu::Atomic::Add(&kcount[s-1], 1);
                        if(power > tolerance)
                        {
                            Gpu::Atomic::Add(&ps_map[s-1], power);   
                        }

                        break;
                    }

                    // If you're at the largest bin
                    else if (kmag == kiso[N/2])
                    { 
                        Real power = (std::pow(field_ptr(i, j, k, component).real(), 2.0) 
                                    + std::pow(field_ptr(i, j, k, component).imag(), 2.0));
                        
                        Gpu::Atomic::Add(&kcount[N/2], 1);
                        if(power > tolerance)
                        {
                            Gpu::Atomic::Add(&ps_map[N/2], power);
                        }

                        break;
                    }

                    // If you've reached the largest bin but not been captured
                    else if(s > N/2)
                    { 
                        Print() << iv << "\n";
                        Print() << kmag << "\n";
                        Print() << kiso[s] << "," << kiso[s-1] << "\n";
                        Error("TensorExtraction::print_power_spectrum Part of the spectrum isn't captured.");
                    }

                    // If you haven't found the right bin yet
                    else { continue; }
                }
            }
        });
    }

    ParallelAllReduce::Sum(kcount.data(), static_cast<int>(kcount.size()), ParallelContext::CommunicatorSub());
    ParallelAllReduce::Sum(ps_map.data(), static_cast<int>(ps_map.size()), ParallelContext::CommunicatorSub());

    // Print the power spectrum to a new file in data/
    for(int s=0; s<=N/2; s++)
    {
        power_spec_file.write_data_line({kiso[s], (double)ps_map[s]/kcount[s]});
    }
}

// Finds statistical moment x of given MultiFab
template <class matter_t>
inline Real TensorExtraction<matter_t>::find_field_moment_x(MultiFab &field, const Vector<Real> mean, 
                                             const int moment, const int component)
{
    Real sum = 0.;
    const Real vol = std::pow(N, 3.);

    // Loop over the box to calculate moment x
    for (MFIter mfi(field); mfi.isValid(); ++mfi) 
    {
        Array4<Real> const& field_ptr = field.array(mfi);
        const Box& bx = mfi.fabbox();

        amrex::ParallelFor(bx, [=, &sum] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            sum += std::pow(field_ptr(i, j, k, component) - mean[component], moment);
        });
    }

    ParallelAllReduce::Sum(sum, ParallelContext::CommunicatorSub());

    // Normalise and return moment x
    if(moment == 2) { return sqrt(sum/vol); }
    else { return sum/vol; }
}

// Calculates and prints requested moments (any between 1 and 4)
template <class matter_t>
inline void TensorExtraction<matter_t>::print_tensor_moment(MultiFab &field, const Vector<std::string> names,  
                                             const Vector<int> &moment_orders, SmallDataIO &file, 
                                             const int is_first_step)
{
    // Trap instance where the user requests too large a moment
    for(const auto moment : moment_orders)
    {
        if(moment > 4) 
        { 
            Error("TensorExtraction::print_tensor_moment Chosen moment order has not been implemented");
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
    Vector<int>::const_iterator stdev_itr = std::find(moment_orders.begin(), moment_orders.end(), 2);
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
    if(is_first_step) { file.write_header_line(headers); }
    file.write_time_data_line(data_to_print);
}

// Extraction routine for hplus and hcross only, called in Example::derive
template <class matter_t>
inline void TensorExtraction<matter_t>::derive(const MultiFab &source, MultiFab &out, int dcomp)
{
    BL_PROFILE("TensorExtraction::derive");

    // Extract MultiFab ingredients from state
    BoxArray sba = source.boxArray();
    DistributionMapping sdm = source.DistributionMap();
    MultiFab hij_x(sba, sdm, 6, 0);

    // Copy the spatial metric from the state
    Copy(hij_x, source, c_h11, lut[0][0], 1, 0);
    Copy(hij_x, source, c_h12, lut[0][1], 1, 0);
    Copy(hij_x, source, c_h13, lut[0][2], 1, 0);
    Copy(hij_x, source, c_h22, lut[1][1], 1, 0);
    Copy(hij_x, source, c_h23, lut[1][2], 1, 0);
    Copy(hij_x, source, c_h33, lut[2][2], 1, 0);

    // Undo the normalisation and BSSN-CPT conversion
    for (int l=0; l<3; l++) { hij_x.plus(-1., lut[l][l], 1); }
    hij_x.mult(1./norm);

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
    cMultiFab hij_k(kba, kdm, 6, 0);

    // Set up the FFT
    IntVect x_domain_high(N-1, N-1, N-1);
    Box x_domain(domain_low, x_domain_high);
    FFT::R2C<Real> tensor_fft(x_domain, FFT::Info().setBatchSize(hij_k.nComp()));

    // Perform the fft
    tensor_fft.forward(hij_x, hij_k);

    // Normalise the fft (fftw style)
    for(int comp = 0; comp < 6; comp++)
    {
        hij_k.mult(std::pow(N, -3.), comp, 1); 
    }

    // Loop to extract the Fourier-space mode functions
    for (MFIter mfi(hij_k); mfi.isValid(); ++mfi) 
    {
        const Box& bx = mfi.fabbox();

        // Make a pointer to the mode functions at this MF box
        Array4<GpuComplex<Real>> const& hs_ptr = hs_k.array(mfi);
        Array4<GpuComplex<Real>> const& hij_ptr = hij_k.array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            IntVect iv{i, j, k};
            Vector<Real> mhat(3, 0.);
            Vector<Real> nhat(3, 0.);

            mhat = calculate_basis_vector(iv, 0);
            nhat = calculate_basis_vector(iv, 1);

            Real eplus = 0.;
            Real ecross = 0.;

            // Find basis tensors and do the Fourier trick
            for (int l=0; l<3; l++) for (int p=0; p<3; p++)
            {
                eplus = mhat[l]*mhat[p] - nhat[l]*nhat[p];
                ecross = mhat[l]*nhat[p] + nhat[l]*mhat[p];

                hs_ptr(i, j, k, 0) += (hij_ptr(i, j, k, lut[l][p]) * eplus)/std::sqrt(2.);
                hs_ptr(i, j, k, 1) += (hij_ptr(i, j, k, lut[l][p]) * ecross)/std::sqrt(2.);
            }
        });
    }

    apply_nyquist_conditions(hs_k);

    // Make a multifab to store config space mode functions
    // Need to use out to make these ingredients??
    BoxArray xba = out.boxArray();//(x_domain); //
    DistributionMapping xdm = out.DistributionMap();//(xba); //
    MultiFab hs_x(xba, xdm, 2, 0);

    // Fourier transform
    FFT::R2C<Real> mode_function_fft(x_domain, FFT::Info().setBatchSize(hs_k.nComp()));
    mode_function_fft.backward(hs_k, hs_x);

    // Apply physical normalisation
    hs_x.mult(norm);

    for (MFIter mfi(hs_x); mfi.isValid(); ++mfi) 
    {
        Array4<Real> const& out_ptr = out.array(mfi);
        Array4<Real> const& hx_ptr = hs_x.array(mfi);

        const Box& bx = mfi.fabbox();
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            const IntVect iv{i, j, k};
            bool in_ghost_index = is_ghost_index(iv);
            if(!in_ghost_index)
            {
                out_ptr(iv, dcomp) = hx_ptr(i, j, k, 0);
                out_ptr(iv, dcomp + 1) = hx_ptr(i, j, k, 1);
            }
        });
    }
}

// Extraction routine called in specificPostTimeStep
template <class matter_t>
inline void TensorExtraction<matter_t>::extract_from_Weyl4(const MultiFab &state, const std::string data_path, const Real dt,  
                                 const Real cur_time, const int restart_time, const int first_step, const int plot_int)
{
    BL_PROFILE("TensorExtraction::extract_from_Weyl4");

    // Extract MultiFab ingredients from state
    BoxArray sba = state.boxArray();
    DistributionMapping sdm = state.DistributionMap();
    MultiFab B_x(sba, sdm, 6, 0);

    const auto &state_arrays = state.const_arrays();
    const auto &B_arrays = B_x.arrays();

    amrex::ParallelFor(B_x, 
        [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k) noexcept
        {
            const auto vars = load_vars<Vars>(state_arrays[box_no].cellData(i, j, k));
            const auto d1   = m_deriv.template diff1<Vars>(i, j, k, state_arrays[box_no]);
            const auto d2   = m_deriv.template diff2<Diff2Vars>(i, j, k, state_arrays[box_no]);

            // Compute the inverse metric
            const auto h_UU  = TensorAlgebra::compute_inverse_sym(vars.h);
            const auto chris = TensorAlgebra::compute_christoffel(d1.h, h_UU);

            // Compute the spatial volume element
            const auto epsilon3_LUU = compute_epsilon3_LUU(vars, h_UU);

            // Compute the E and B fields
            EBFields_t<Real> ebfields =
                compute_EB_fields(vars, d1, d2, epsilon3_LUU, h_UU, chris);

            // Add in matter terms to E and B fields
            add_matter_EB(ebfields, vars, d1, epsilon3_LUU, h_UU, chris);

            for(int m=0; m<3; m++) for(int n=m; n<3; n++)
            {
                B_arrays[box_no](i, j, k, lut[m][n]) = ebfields.B[m][n];
            }
        });

    // Set up the problem domain in Fourier space
    // And impose that MPI ranks only slice along the i index (for Nyquist conditions)
    IntVect domain_low(0, 0, 0);
    IntVect k_domain_high(N/2, N-1, N-1);
    Box k_domain(domain_low, k_domain_high);
    Array< bool, AMREX_SPACEDIM > const &slicing{true, false, false};
    BoxArray kba = decompose(k_domain, ParallelContext::NProcsAll(), slicing);
    DistributionMapping kdm(kba);

    // Set up the arrays to store the Fourier data sets
    cMultiFab B_k(kba, kdm, 6, 0);
    cMultiFab hdot_k(kba, kdm, 2, 0);

    // Set up the FFT
    IntVect x_domain_high(N-1, N-1, N-1);
    Box x_domain(domain_low, x_domain_high);
    FFT::R2C<Real> tensor_fft(x_domain, FFT::Info().setBatchSize(B_k.nComp()));

    // Perform the fft
    tensor_fft.forward(B_x, B_k);

    // Normalise the fft (fftw style)
    for(int comp = 0; comp < 6; comp++)
    {
        B_k.mult(std::pow(N, -3.), comp, 1); 
    }

    std::string Filename = "/nfs/st01/hpc-gr-epss/eaf49/GRTeclyn-dump/hs-k-extr";
    int time_step = cur_time/dt;

    const auto &Bk_arrays = B_k.arrays();
    const auto &hdotk_arrays = hdot_k.arrays();

    // Loop to extract the Fourier-space mode functions
    amrex::ParallelFor(B_k, 
    [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k) noexcept
    {
        // Turn B -> h'ij by inverting the curl
        const auto vars = load_vars<Vars>(state_arrays[box_no].cellData(i, j, k));
        const auto epsilon3_LLL = compute_epsilon3_LLL(vars);
        
        IntVect iv{i, j, k};
        Vector<Real> k_vec(iv.size());
        for(int i=0; i<iv.size(); i++) { k_vec[i] = static_cast<Real>(iv[i]) * 2. * M_PI / L; }
        Real kmag_sq = (i*i + j*j + k*k) * 2 * M_PI / L;

        Vector<Real> mhat(3, 0.);
        Vector<Real> nhat(3, 0.);

        mhat = calculate_basis_vector(iv, 0);
        nhat = calculate_basis_vector(iv, 1);

        Real eplus = 0.;
        Real ecross = 0.;

        // Find basis tensors and do the Fourier trick
        for (int l=0; l<3; l++) for (int p=0; p<3; p++)
        {
            for(int m=0; m<3; m++)
            {
                Bk_arrays[box_no](i, j, k, lut[l][p]) *= k_vec[m] * epsilon3_LLL[l][p][m] / kmag_sq;
                Bk_arrays[box_no](i, j, k, lut[l][p]) = swap_real_imag_parts(Bk_arrays[box_no](i, j, k, lut[l][p]));
            }

            eplus = mhat[l]*mhat[p] - nhat[l]*nhat[p];
            ecross = mhat[l]*nhat[p] + nhat[l]*mhat[p];

            hdotk_arrays[box_no](i, j, k, 0) += (Bk_arrays[box_no](i, j, k, lut[l][p]) * eplus)/std::sqrt(2.);
            hdotk_arrays[box_no](i, j, k, 1) += (Bk_arrays[box_no](i, j, k, lut[l][p]) * ecross)/std::sqrt(2.);
        }
    });

    apply_nyquist_conditions(hdot_k);

    // Find the binned PS for each mode function and print to data/
    if((m_params.calc_binned_power_spectrum) && (time_step % plot_int == 0)) 
    {
        std::string spec_path = make_subdirectory(data_path, "spectra", first_step);
        Vector<std::string> filenames(2, "");

        for(int comp = 0; comp < hdot_k.nComp(); comp++)
        {
            filenames[comp] = spec_path+"Weyl-spectrum-comp-"+std::to_string(comp)+"-time-";
            SmallDataIO spectrum_file(filenames[comp], dt, cur_time, restart_time, SmallDataIO::NEW, first_step, ".dat");
            print_power_spectrum(hdot_k, spectrum_file, comp);
        }
    }
}

// Extraction routine called in specificPostTimeStep
template <class matter_t>
inline void TensorExtraction<matter_t>::extract(const MultiFab &state, const std::string data_path, const Real dt,  
                                 const Real cur_time, const int restart_time, const int first_step, const int plot_int)
{
    BL_PROFILE("TensorExtraction::extract");

    // Extract MultiFab ingredients from state
    BoxArray sba = state.boxArray();
    DistributionMapping sdm = state.DistributionMap();
    MultiFab Aij_x(sba, sdm, 6, 0);

    // Copy the spatial metric from the state
    Copy(Aij_x, state, c_A11, lut[0][0], 1, 0);
    Copy(Aij_x, state, c_A12, lut[0][1], 1, 0);
    Copy(Aij_x, state, c_A13, lut[0][2], 1, 0);
    Copy(Aij_x, state, c_A22, lut[1][1], 1, 0);
    Copy(Aij_x, state, c_A23, lut[1][2], 1, 0);
    Copy(Aij_x, state, c_A33, lut[2][2], 1, 0);

    // Undo the normalisation and BSSN-CPT conversion
    //for (int l=0; l<3; l++) { hij_x.plus(-1., lut[l][l], 1); }
    Aij_x.mult(2./norm);

    // Set up the problem domain in Fourier space
    // And impose that MPI ranks only slice along the i index (for Nyquist conditions)
    IntVect domain_low(0, 0, 0);
    IntVect k_domain_high(N/2, N-1, N-1);
    Box k_domain(domain_low, k_domain_high);
    Array< bool, AMREX_SPACEDIM > const &slicing{true, false, false};
    BoxArray kba = decompose(k_domain, ParallelContext::NProcsAll(), slicing);
    DistributionMapping kdm(kba);

    // Set up the arrays to store the Fourier data sets
    cMultiFab As_k(kba, kdm, 2, 0);
    cMultiFab Aij_k(kba, kdm, 6, 0);

    // Set up the FFT
    IntVect x_domain_high(N-1, N-1, N-1);
    Box x_domain(domain_low, x_domain_high);
    FFT::R2C<Real> tensor_fft(x_domain, FFT::Info().setBatchSize(Aij_k.nComp()));

    // Perform the fft
    tensor_fft.forward(Aij_x, Aij_k);

    // Normalise the fft (fftw style)
    for(int comp = 0; comp < 6; comp++)
    {
        Aij_k.mult(std::pow(N, -3.), comp, 1); 
    }

    std::string Filename = "/nfs/st01/hpc-gr-epss/eaf49/GRTeclyn-dump/As-k-extr";
    int time_step = cur_time/dt;

    // Loop to extract the Fourier-space mode functions
    for (MFIter mfi(Aij_k); mfi.isValid(); ++mfi) 
    {
        const Box& bx = mfi.fabbox();

        // Make a pointer to the mode functions at this MF box
        Array4<GpuComplex<Real>> const& As_ptr = As_k.array(mfi);
        Array4<GpuComplex<Real>> const& Aij_ptr = Aij_k.array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            IntVect iv{i, j, k};
            Vector<Real> mhat(3, 0.);
            Vector<Real> nhat(3, 0.);

            mhat = calculate_basis_vector(iv, 0);
            nhat = calculate_basis_vector(iv, 1);

            Real eplus = 0.;
            Real ecross = 0.;

            // Find basis tensors and do the Fourier trick
            for (int l=0; l<3; l++) for (int p=0; p<3; p++)
            {
                eplus = mhat[l]*mhat[p] - nhat[l]*nhat[p];
                ecross = mhat[l]*nhat[p] + nhat[l]*mhat[p];

                As_ptr(i, j, k, 0) += (Aij_ptr(i, j, k, lut[l][p]) * eplus)/std::sqrt(2.);
                As_ptr(i, j, k, 1) += (Aij_ptr(i, j, k, lut[l][p]) * ecross)/std::sqrt(2.);
            }
        });
    }

    apply_nyquist_conditions(As_k);

    // Find the binned PS for each mode function and print to data/
    if((m_params.calc_binned_power_spectrum) && (time_step % plot_int == 0)) 
    {
        std::string spec_path = make_subdirectory(data_path, "spectra", first_step);
        Vector<std::string> filenames(2, "");

        for(int comp = 0; comp < As_k.nComp(); comp++)
        {
            filenames[comp] = spec_path+"TT-spectrum-comp-"+std::to_string(comp)+"-time-";
            SmallDataIO spectrum_file(filenames[comp], dt, cur_time, restart_time, SmallDataIO::NEW, first_step, ".dat");
            print_power_spectrum(As_k, spectrum_file, comp);
        }
    }

    // Find mode functions in configuration space if requested
    if(m_params.calc_higher_order_statistics)
    {
        // Make a multifab to store config space mode functions
        BoxArray xba(x_domain);
        DistributionMapping xdm(xba);
        MultiFab As_x(xba, xdm, 2, 0);

        // Fourier transform
        FFT::R2C<Real> mode_function_fft(x_domain, FFT::Info().setBatchSize(As_k.nComp()));
        mode_function_fft.backward(As_k, As_x);

        // Apply physical normalisation
        As_x.mult(norm);

        // Print mode functions if requested
        /*std::string mf_path = make_subdirectory(data_path, "mode-functions", first_step);
        std::string filename = mf_path+"mode-function-"+std::to_string(cur_time/dt);

        for (MFIter mfi(As_x); mfi.isValid(); ++mfi) 
        {
            Array4<Real> const& hx_ptr = As_x.array(mfi);
            const Box& bx = mfi.fabbox();

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                AllPrintToFile(filename) << i + N*(j + N*k) << ",";
                for(int c=0; c<2; c++)
                {
                    AllPrintToFile(filename).SetPrecision(14) << hx_ptr(i, j, k, c) << ",";
                }
                AllPrintToFile(filename) << "\n";
            });
        }*/

        // Calculate and print field moments if requested
        if (m_params.calc_higher_order_statistics)
        {
            SmallDataIO stats_file(data_path+"field-statistics", dt, cur_time, 
                                    restart_time, SmallDataIO::APPEND, first_step, ".dat");

            if (!m_params.orders.empty())
            {
                Vector<std::string> names{"hplus","hcross"};
                print_tensor_moment(As_x, names, m_params.orders, stats_file, first_step);
            }
        }
    }
}

#endif
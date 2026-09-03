/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef DERIVEDVARIABLES_IMPL_HPP_
#define DERIVEDVARIABLES_IMPL_HPP_

void DerivedVariables::set_up(int a_state_index)
{
    int num_ghosts = 0;

    auto &derive_lst     = amrex::AmrLevel::get_derive_lst();
    const auto &desc_lst = amrex::AmrLevel::get_desc_lst();

    const auto &field_names = DerivedVariables::var_names;

    // Add Constraints to the derive list
    derive_lst.add(
        DerivedVariables::name, amrex::IndexType::TheCellType(),
        static_cast<int>(field_names.size()), field_names,
        DerivedVariables::compute_mf, [=](const amrex::Box &box)
        { return amrex::grow(box, num_ghosts); }, &amrex::cell_quartic_interp);

    derive_lst.addComponent(DerivedVariables::name, desc_lst, a_state_index, 0,
                            NUM_VARS);
}

/* Main functions */

// Extract R and hs in configuration space from the BSSN variables
inline void DerivedVariables::extract_hs_and_R(amrex::MultiFab &hs,
                                               amrex::MultiFab &R,
                                               const amrex::MultiFab &state)
{
    // Extract amrex::MultiFab ingredients from state
    amrex::BoxArray sba            = state.boxArray();
    amrex::DistributionMapping sdm = state.DistributionMap();
    if (sba != hs.boxArray() || sdm != hs.DistributionMap())
    {
        amrex::Error("DerivedVariables::extract_hs_and_R "
                     "source and output BA or SDM do not match");
    }

    // 0: scalar field
    // 1: conformal factor
    amrex::MultiFab scalars_x(sba, sdm, 2, 0);
    amrex::MultiFab gij_x(sba, sdm, 6, 0);

    // Copy the spatial metric from the state
    Copy(gij_x, state, c_h11, InflatonUtils::look_up_table[0][0], 1, 0);
    Copy(gij_x, state, c_h12, InflatonUtils::look_up_table[0][1], 1, 0);
    Copy(gij_x, state, c_h13, InflatonUtils::look_up_table[0][2], 1, 0);
    Copy(gij_x, state, c_h22, InflatonUtils::look_up_table[1][1], 1, 0);
    Copy(gij_x, state, c_h23, InflatonUtils::look_up_table[1][2], 1, 0);
    Copy(gij_x, state, c_h33, InflatonUtils::look_up_table[2][2], 1, 0);

    constexpr int phi_component = 0;
    constexpr int chi_component = 1;
    Copy(scalars_x, state, c_phi, phi_component, 1, 0);
    Copy(scalars_x, state, c_chi, chi_component, 1, 0);

    // Find background quantities needed to extract \cal R
    const int vol               = std::pow(params().N, 3);
    const amrex::Real K_bar     = state.sum(c_K) / vol;
    const amrex::Real alpha_bar = state.sum(c_lapse) / vol;
    const amrex::Real Pi_bar    = state.sum(c_Pi) / vol;
    const amrex::Real phi_bar   = state.sum(c_phi) / vol;
    const amrex::Real chi_bar   = state.sum(c_chi) / vol;

    // Remove background from scalar field
    scalars_x.plus(-phi_bar, phi_component, 1);
    scalars_x.plus(-chi_bar, chi_component, 1);
    scalars_x.mult(1. / m_utils.calculate_norm());

    // Undo the normalisation and BSSN-CPT conversion
    for (int l = 0; l < 3; l++)
    {
        gij_x.plus(-1., InflatonUtils::look_up_table[l][l], 1);
    }
    gij_x.mult(1. / m_utils.calculate_norm());

    // Set up the problem domain in Fourier space
    // And impose that MPI ranks only slice along the i index (for Nyquist
    // conditions)
    amrex::IntVect domain_low(0, 0, 0);
    amrex::IntVect k_domain_high(params().N / 2, params().N - 1,
                                 params().N - 1);
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
    amrex::IntVect x_domain_high(params().N - 1, params().N - 1,
                                 params().N - 1);
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
    for (int comp = 0; comp < 6; comp++)
    {
        gij_k.mult(std::pow(params().N, -3.), comp, 1);
    }
    for (int comp = 0; comp < 2; comp++)
    {
        scalars_k.mult(std::pow(params().N, -3.), comp, 1);
    }

    const auto &hs_arrs      = hs_k.arrays();
    const auto &gij_arrs     = gij_k.arrays();
    const auto &scalars_arrs = scalars_k.arrays();
    const auto &R_k_arrs     = R_k.arrays();

    // Local copy so the kernel captures config by value, not via the host
    // `this` pointer
    const InflatonUtils cfg = m_utils;

    amrex::ParallelFor(
        gij_k,
        [=] AMREX_GPU_DEVICE(int bx, int i, int j, int k)
        {
            amrex::IntVect iv{i, j, k};
            amrex::Real kmag = cfg.get_kmag(iv);

            if (iv != amrex::IntVect{0, 0, 0})
            {
                const auto [eplus, ecross] =
                    cfg.calculate_polarisation_tensors(iv);

                // Find basis tensors and do the Fourier trick
                for (int l = 0; l < 3; l++)
                    for (int p = 0; p < 3; p++)
                    {
                        hs_arrs[bx](i, j, k, 0) +=
                            (gij_arrs[bx](i, j, k,
                                          InflatonUtils::look_up_table[l][p]) *
                             eplus(l, p)) /
                            2.;
                        hs_arrs[bx](i, j, k, 1) +=
                            (gij_arrs[bx](i, j, k,
                                          InflatonUtils::look_up_table[l][p]) *
                             ecross(l, p)) /
                            2.;
                    }

                // Calculate the TT and scalar-(vector) components of the
                // metric, by reconstructing hij and subtracting it from
                // \tilde{gamma}_ij
                Tensor::Rank2 hij_re, hSV_im, hij_im, hSV_re;
                for (int l = 0; l < 3; l++)
                    for (int p = 0; p < 3; p++)
                    {
                        hij_re(l, p) =
                            (hs_arrs[bx](i, j, k, 0).real() * eplus(l, p) +
                             hs_arrs[bx](i, j, k, 1).real() * ecross(l, p));
                        hij_im(l, p) =
                            (hs_arrs[bx](i, j, k, 0).imag() * eplus(l, p) +
                             hs_arrs[bx](i, j, k, 1).imag() * ecross(l, p));

                        hSV_re(l, p) =
                            gij_arrs[bx](i, j, k,
                                         InflatonUtils::look_up_table[l][p])
                                .real() -
                            hij_re(l, p);
                        hSV_im(l, p) =
                            gij_arrs[bx](i, j, k,
                                         InflatonUtils::look_up_table[l][p])
                                .imag() -
                            hij_im(l, p);
                    }

                // Extract R according to the scheme detailed in
                // Appendix B (Eq. B1) of arxiv:2502.06783, using hSV as the
                // spatial metric instead of \tilde{gamma}_ij
                if (cfg.m_params.scalar_init)
                {
                    // Find the unitful k vector
                    amrex::GpuArray<amrex::Real, 3> iv_k{
                        static_cast<amrex::Real>(iv[0]),
                        static_cast<amrex::Real>(iv[1]),
                        static_cast<amrex::Real>(iv[2])};
                    iv_k[1] = cfg.invert_index_with_sign(iv_k[1]);
                    iv_k[2] = cfg.invert_index_with_sign(iv_k[2]);

                    for (auto &k_comp : iv_k)
                    {
                        k_comp *= 2. * amrex::Math::pi<amrex::Real>() /
                                  cfg.m_params.L;
                    }
                    amrex::GpuComplex<amrex::Real> Phi = 0;

                    // converstion from chi and gamma_ij -> Phi
                    for (int l = 0; l < 3; l++)
                        for (int p = 0; p < 3; p++)
                        {
                            Phi += amrex::GpuComplex<amrex::Real>{
                                (iv_k[l] * iv_k[p] * hSV_re(l, p)) /
                                    std::pow(kmag, 2.),
                                (iv_k[l] * iv_k[p] * hSV_im(l, p)) /
                                    std::pow(kmag, 2.)};
                        }
                    Phi *= 1. / 4.;
                    Phi += 0.5 * (scalars_arrs[bx](i, j, k, chi_component));

                    // Combine the above to find R(k)
                    R_k_arrs[bx](i, j, k, 0) =
                        Phi - ((K_bar / 3.) *
                               scalars_arrs[bx](i, j, k, phi_component) /
                               alpha_bar / Pi_bar);
                }
            }
        });

    amrex::Gpu::streamSynchronize();

    // Prepare to IFT the polarisation fields and R field
    m_utils.apply_nyquist_conditions(hs_k);
    m_utils.apply_nyquist_conditions(R_k);

    // Fourier transform
    mode_fn_fft.backward(hs_k, hs);
    R_fft.backward(R_k, R);

    // Apply physical normalisation
    hs.mult(m_utils.calculate_norm());
    R.mult(m_utils.calculate_norm());
}

// Put R and hs into plotfiles
// DeriveFuncMF callback: build an extractor and fill the plotfile output with
// R, hplus, hcross. src_mf arrives already FillPatch-ed by the framework.
inline void DerivedVariables::compute_mf(amrex::MultiFab &out_mf, int dcomp,
                                         int /*ncomp*/,
                                         const amrex::MultiFab &src_mf,
                                         const amrex::Geometry & /*geomdata*/,
                                         amrex::Real /*time*/,
                                         const int * /*bcrec*/, int /*level*/)
{
    BL_PROFILE("DerivedVariables::compute_mf");

    // Make a multifab to store config space mode functions
    amrex::BoxArray oba            = out_mf.boxArray();
    amrex::DistributionMapping odm = out_mf.DistributionMap();
    amrex::MultiFab hs_x(oba, odm, 2, 0);
    amrex::MultiFab R_x(oba, odm, 1, 0);
    hs_x.setVal(0.0);
    R_x.setVal(0.0);

    DerivedVariables extractor;
    extractor.extract_hs_and_R(hs_x, R_x, src_mf);

    const auto &hs_arrs  = hs_x.arrays();
    const auto &R_arrs   = R_x.arrays();
    const auto &out_arrs = out_mf.arrays();

    amrex::ParallelFor(hs_x,
                       [=] AMREX_GPU_DEVICE(int bx, int i, int j, int k)
                       {
                           const amrex::IntVect iv{i, j, k};
                           out_arrs[bx](iv, dcomp) = R_arrs[bx](i, j, k);
                           out_arrs[bx](iv, dcomp + 1) =
                               hs_arrs[bx](i, j, k, 0);
                           out_arrs[bx](iv, dcomp + 2) =
                               hs_arrs[bx](i, j, k, 1);
                       });

    amrex::Gpu::streamSynchronize();
}

#endif /* DERIVEDVARIABLES_IMPL_HPP_ */

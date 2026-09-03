# The inflaton field example

In the following section, we describe how to run the inflaton field example. This example implements the classical stochastic initial conditions used in lattice cosmology simulations of inflation. This example produces Gaussian random initial conditions for scalar and gravitational-wave fluctuations which satisfy the constraints at first order, on top of the chosen inflationary background. The user then has the option to calculate the scalar perturbations, in the form of the gauge-invariant curvature perturbation $\mathcal{R}$, and the gravitational waves, in the form of the plus and cross polarisation fields, as diagnostic quantities, and to print these fields to plot files.

__(!!) As this example relies on the use of a Fast Fourier Transform (FFT), the user must be particularly careful about their use of AMR when running this example (!!)__

## Physical scenario

This page describes how to run the inflaton field example using a parameter file similar to [the default parameter file](https://github.com/the-florist/GRTeclyn-workspace/blob/example/InflatonField/Examples/InflatonField/params.txt) (UPDATEME) given in the base directory. 

First, a few qualifiers on the use of this example...

### Initial data 

The initialisation for this code comes in two parts. The user must supply a background inflationary model, determined by the initial inflaton field mean, velocity and mass parameter. The background or ``average`` values of the BSSN variables are then determined from the appropriate cosmological quantities, via the BSSN-CPT correspondence. This correspondence is derived for scalar and tensor perturbations in [arXiv:2401.08530](https://arxiv.org/abs/2401.08530) and [arXiv:2412.19731](https://arxiv.org/abs/2412.19731). For the background, this correspondence is given by 

$$ 
\chi = \frac{1}{a^2},\ \tilde{\gamma}_{ij} = \delta_{ij},\ K = - 3 H,\ A_{ij} = 0,
$$

where $a$ is the scale factor and $H$ is the Hubble parameter.

Once the background is initialised, perturbation fields are constructed for the BSSN variables in Fourier space. These perturbation fields are then inverse Fourier-transformed and added to the appropriate background BSSN quantites. The perturbation fields are Gaussian random and obey the analytic solution to the first-order equation of motion for scalar and tensor perturbations, respectively, during inflation, and satisfy the constraints at first order.
The scalar field and the two gravitational-wave polarisation fields both satisfy the Mukhanov-Sasaki equation at first order:

$$
f'' + \left(k^2 - \frac{a''}{a}\right)f = 0.
$$

The analytic solution to this equation, which assumes a Bunch-Davies vacuum state in the far past, gives the following polarisation field mode function at the start of the simulation:

$$
h(k) = \frac{1}{M_p a_0 \sqrt{2k}}\left(1 - i\frac{a_0 H_0}{k}\right) e^{-ik/(a_0 H_0)}
$$

where $a_0$ and $H_0$ refer to the initial scale factor and Hubble parameter values. The initial mode function for $\mathcal{R}$ is simply 

$$
\mathcal{R}(k) = \frac{1}{\sqrt{2 \epsilon_1}} f(k)
$$

where $\epsilon_1$ is the first slow-roll parameter. See [arXiv:2401.08530](https://arxiv.org/abs/2401.08530) and [arXiv:2412.19731](https://arxiv.org/abs/2412.19731) for a full derivation of these relations, along with expressions for the velocity mode function solutions.

We first break these two analytic solutions into their modulus and phase, and then initialise the perturbation fields as 

$$
f(\mathbf{k}) = |f(k)|\ \sqrt{-2 \ln u_1}\ e^{i \left[ \varphi_f(k) + 2 \pi u_2 \right]},
$$

where $|f(k)|$ and $\varphi_f(k)$ are the modulus and phase of the analytic mode function, and $u_1$ and $u_2$ are two random numbers drawn uniformly from the open interval $(0, 1)$. It can be shown that a field constructed in this manner in Fourier space will be Gaussian random in configuration space. 

Rather than drawing $u_1$ and $u_2$ from a sequential random number generator, we obtain them by hashing a key built from `init.random_seed`, the current grid location. This ensures that for the same random seed, the same set of draws are made no matter how the domain is divided between MPI ranks or threads. Note also that the amplitude and velocity mode functions for each field share a single pair of draws, since they describe the same stochastic realisation and its time derivative, rather than two independent ones.

The user can then choose to apply a window function to the initial power spectrum, which can be useful in reducing high-frequency noise while retaining a continuous spectrum at low modes. The window function is applied where `init.use_window = 1`, and is given by

$$
W(k) = \frac{1}{2} \left[ 1 - \tanh \left( \frac{L}{\Delta} \left( k - k_s \right) \right) \right], \qquad k_s = \frac{\sqrt{3}\, \pi N_w}{10 L},
$$

where the window function width is set by $\Delta$, and wherre $N_w$ is `N_coarse` if it has been supplied, and the grid resolution $N$ otherwise. This rolls the spectrum smoothly off above $k_s$, removing power close to the grid scale which would otherwise be poorly resolved. Because $k_s$ can be tied to the coarsest resolution of interest rather than to $N$, the same physical initial data can be represented on grids of different resolution, which is what makes convergence testing possible.

The two sectors are then assembled from these random fields. In the tensor sector, two independent polarisation amplitudes are contracted with the plus and cross polarisation tensors, which are built from an orthonormal pair of basis vectors transverse to $\mathbf{k}$, to give the perturbations to $\tilde{\gamma}_{ij}$ and $A_{ij}$. The parameter `init.alpha` rotates this pair of basis vectors within the transverse plane. In the scalar sector, the single realisation of $\mathcal{R}$ and its derivative is converted into perturbations of $\phi$, $\Pi$, $\chi$ and $K$ which satisfy the constraints at first order. Both constructions follow from the BSSN-CPT correspondence, and we again refer the reader to [arXiv:2401.08530](https://arxiv.org/abs/2401.08530) and [arXiv:2412.19731](https://arxiv.org/abs/2412.19731) for their derivation.

Finally, the fields are transformed back to configuration space using AMReX's real-to-complex transform, `amrex::FFT::R2C`, and added to the background. AMReX's transforms are unnormalised, so we apply the normalisation

$$
\left( \frac{\sqrt{2 \pi}}{L} \right)^3
$$

to each field after the backward transform. This is the convention in which the mode functions quoted above are the physical ones, and the same factor is applied in reverse by the diagnostic extraction described in the next section, so that a field which is written out and read back in returns the spectrum it started with. It should be changed with care, as the two must agree.

### Derived variables

Alongside the evolution variables, this example can output three diagnostic fields: the gauge-invariant curvature perturbation $\mathcal{R}$, and the two gravitational-wave polarisation fields, printed as `hplus` and `hcross`. These are AMReX derived variables, which are reconstructed on demand whenever a plot file is written, and they are requested by adding `InflationFields` to `amr.derive_plot_vars`. Note that because they are reconstructed using a Fourier transform over the whole domain, they cannot be used together with AMR, and the code will stop with an error if these variables are requested while `amr.max_level` is non-zero.

The task at extraction is essentially the reverse of the initialisation: given only the evolved BSSN variables, we wish to separate the scalar and tensor degrees of freedom at first order. We first remove the background from the conformal metric,

$$
\delta\tilde{\gamma}_{ij} = \tilde{\gamma}_{ij} - \delta_{ij},
$$

and likewise from the scalar field and the conformal factor, giving $\delta\phi$ and $\delta\chi$. The background values $\bar{\phi}$ and $\bar{\chi}$, along with $\bar{K}$, $\bar{\alpha}$ and $\bar{\Pi}$, are computed as averages over the grid. These fields are then transformed into Fourier space.

The transverse-traceless part of the metric perturbation is obtained by decomposing $\delta\tilde{\gamma}_{ij}$ onto the same polarisation basis used to construct the initial data,

$$
h^{TT}_{ij} = \sum_s \text{e}^s_{ij} \delta\tilde{\gamma}_s, \qquad \delta\tilde{\gamma}_s = \frac{1}{2}\sum_{i,j} \delta\tilde{\gamma}_{ij}\, \text{e}_s^{ij},
$$

where the factor of one half in the projection follows from the normalisation of the polarisation tensors. The two amplitudes $\delta\tilde{\gamma}_+$ and $\delta\tilde{\gamma}_{\times}$ are the fields printed as `hplus` and `hcross`. Since the polarisation tensors are themselves transverse and trace-free, $h^{TT}_{ij}$ is transverse and trace-free by construction. The remaining scalar and vector part of the metric perturbation then follows as

$$
h^{SV}_{ij} = \delta\tilde{\gamma}_{ij} - h^{TT}_{ij}.
$$

The curvature perturbation can then be reconstructed from the BSSN scalars together with $h^{SV}_{ij}$, following the scheme derived in Appendix B of [arXiv:2502.06783](https://arxiv.org/abs/2502.06783):

$$
\mathcal{R}(\textbf{k}) = \frac{1}{2}\delta\chi(\textbf{k}) - \frac{\bar{K}}{3 \bar{\alpha} \bar{\Pi}}\delta\phi(\textbf{k}) + \frac{1}{4}\sum_{i,j}\frac{k^i k^j}{k^2}h^{SV}_{ij}(\textbf{k}),
$$

where $\Pi = \partial_t \phi$. The three diagnostic fields are then transformed back into configuration space, with the same Fourier normalisation convention that was applied to the initial data, so that the spectra measured here can be compared directly against those used to initialise the run.


## Computational set-up

The example is compiled in the same way as any other GRTeclyn example, as described in [Building on CPUs](building_cpus.md) and [Building and Running on GPUs](building_gpus.md). Navigate to the example directory and build:

```bash
cd GRTeclyn/Examples/InflatonField
make -j 4
```

This example differs from the others in that it requires a Fast Fourier Transform library. `USE_FFT = TRUE` is already set in the example's `GNUmakefile`, but FFTW must be available on your system or loaded as a module file. The resulting executable takes the parameter file as its first argument:

```bash
mpiexec -n 4 ./InflatonField3d.gnu.MPI.ex params.txt
```

There are two groups of parameters which must be changed with particular care in this example.
The first is `amr.max_level`, which should be left at zero if the example-specific diagnostics are being used. AMR could in principle be used with this example if there is only one level on the initial slice, however this use for this example has never been tested before. For the same reason the grid must be cubic and periodic, so `amr.n_cell` and `geometry.prob_extent` must be equal in all three directions and `geometry.is_periodic` must be set in all three directions. These conditions are checked at start-up.

The second is the gauge. The initial data is by necessity (?) constructed in generalised synchronous gauge, in which the lapse is unity and the shift vanishes. In the default example the evolution preserves this gauge choice, and the gauge driver is switched off entirely by setting

```
gauge.lapse_coeff = 0.0
gauge.lapse_advec_coeff = 0.0
gauge.shift_Gamma_coeff = 0.0
gauge.shift_advec_coeff = 0.0
gauge.eta = 0.0
```

### Convergence testing

Because the initial conditions in this example are stochastic, some care is needed when running a convergence test. The perturbation fields are generated mode by mode on the Fourier grid, so simply changing $N$ changes both the sampling of the box and the set of modes which the grid can represent. Two runs at different resolutions would then start from genuinely different realisations, and any difference between them would reflect the change of realisation rather than the resolution. The parameters `N_fine` and `N_coarse` exist to remove both of these effects, and both should be set whenever a suite of convergence test runs is to be compared.

`N_fine` fixes the grid on which the realisation is generated. When it is larger than $N$, the random field is constructed at resolution `N_fine` and then sampled down onto the simulation grid, so that every run in the suite is initialised from the same underlying random field regardless of its own resolution. It should therefore be set to the finest resolution in the suite, and must be at least $N$.
`N_coarse` fixes the cut-off scale $k_s$ of the window function described above. It should be set to the coarsest resolution in the suite, so that the spectrum is cut off at a scale which every run can resolve. In practice we typically choose values of $N$ between 64 and 512.

## Getting extra features

## Guide to parameters

The parameters below are specific to this example. All other parameters used in `params.txt` are shared with the other examples and are described in the [Parameters guide](parameters.md).

* `scalar_field.G_Newton`: The value of G which sets the units of the simulation. In the default case we use Planck units, where $G = 1/\sqrt{8\pi}$.
* `init.background_phi`: The value of the initial inflaton field mean in units of $M_p$.
* `init.background_dphi`: The value of the initial inflaton field derivative, $\Pi = \partial_t\phi$.
* `scalar_field.scalar_mass`: The value of the inflaton mass in units of $M_p$.
* `init.scalar_init`: Whether to add scalar perturbations to the background. Set to `1` to enable and `0` to disable.
* `init.tensor_init`: Whether to add gravitational-wave perturbations to the background. Set to `1` to enable and `0` to disable.
* `init.random_seed`: The seed used to generate the random draws.
* `init.use_window`: Whether to apply the window function to the initial spectrum.
* `init.Delta`: The width of the window function, which is applied as $L/\Delta$.
* `init.alpha`: The angle, in degrees, through which the pair of polarisation basis vectors is rotated within the plane transverse to $\textbf{k}$.
* `N_fine`: The grid resolution on which the random realisation is generated, before being sampled down onto the simulation grid. Defaults to $N$, and must be at least $N$. Set this to the finest resolution in a convergence test suite.
* `N_coarse`: The resolution used to set the window function cut-off scale $k_s$. Defaults to $N$, and must be at most $N$. Set this to the coarsest resolution in a convergence test suite.

To use the derived variables associated with this example, add `InflatonFields` to the parameter `amr.derive_plot_vars`. Also, the `ccz4.min_chi` and `ccz4.min_lapse` parameters should be set to very small non-zero values, as $\chi$ is related to the scale factor and so will naturally shrink exponentially in time. We recommend the user use the BSSN formulation without constraint damping, as the use of constraint damping can alter the evolution of the perturbation fields in unpredictable and potentially unphysical ways.

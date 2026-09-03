# The inflaton field example

In the following section, we describe how to run the inflaton field example. This example implements the classical stochastic initial conditions used in lattice cosmology simulations of inflation. This example produces Gaussian random initial conditions for scalar and gravitational-wave fluctuations which satisfy the constraints at first order, on top of the chosen inflationary background. The user then has the option to calculate the scalar perturbations, in the form of the gauge-invariant curvature perturbation $\mathcal{R}$, and the gravitational waves, in the form of the plus and cross polarisation fields, as diagnostic quantities, and to print these fields to plot files.

__(!!) As this example relies on the use of a Fast Fourier Transform (FFT), the user must be particularly careful about their use of AMR when running this example (!!)__

## Physical scenario

This page describes how to run the inflaton field example using a parameter file similar to [the default parameter file](https://github.com/the-florist/GRTeclyn-workspace/blob/example/InflatonField/Examples/InflatonField/params.txt) (UPDATEME) given in the base directory. 

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
h(k) = \frac{1}{M_p a_0 \sqrt{2k}}\left(1 - \frac{i a_0 H_0}{k}\right) e^{-ik/(a_0 H_0)}
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



## Computational set-up

### Convergence testing

## Checking the outputs

## Qualifiers on the use of this example

## Getting extra features

## Guide to parameters
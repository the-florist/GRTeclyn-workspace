# The inflaton field example

In the following section, we describe how to run the inflaton field example. This example implements the classical stochastic initial conditions used in lattice cosmology simulations of inflation. This example produces Gaussian random initial conditions for scalar and gravitational-wave fluctuations which satisfy the constraints at first order, on top of the chosen inflationary background. The user then has the option to calculate the scalar perturbations, in the form of the gauge-invariant curvature perturbation $\mathcal{R}$, and the gravitational waves, in the form of the plus and cross polarisation fields, as diagnostic quantities, and to print these fields to plot files.

_As this example relies on the use of a Fast Fourier Transform (FFT), the user must be particularly careful about their use of AMR when running this example_
We expand on this point in the following sections.

## Physical scenario

This page describes how to run the inflaton field example using a parameter file similar to [the default parameter file]() given in the base directory. 

The initialisation for this code comes in two parts. The user must supply a background inflationary model, determined by 

### Initial data 

### Derived variables

## Computational set-up

## Checking the outputs

## Getting extra features
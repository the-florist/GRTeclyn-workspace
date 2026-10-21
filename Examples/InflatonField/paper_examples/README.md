# Paper examples: arXiv:2502.06783

Ready-to-run inputs for the three scenarios of
*Bunch-Davies initial conditions and non-perturbative inflationary dynamics
in numerical relativity* (Y. L. Launay, G. I. Rigopoulos, E. P. S. Shellard,
[arXiv:2502.06783](https://arxiv.org/abs/2502.06783)):

| scenario | potential (`potential.type`) | σ | L (1/H) small / paper | ΔN before start | stop time (1/H) |
|---|---|---|---|---|---|
| quadratic | ½m²φ², H = 10⁻⁵ Mp, ε₁ = 0.01 (`quadratic`) | 3 | 33.5 / 64 | 5 | 8 |
| inflection point (ultra-slow-roll) | Λ = 1.86×10⁻⁶, v = 0.198 Mp, φ_i = −1.2 Mp (`usr`) | 4 | 25.1 / 32 | 3 | 10 |
| strong resonance (monodromy) | ½m²φ² + windowed cosine, α = 10, b = 1.32, β = 2 (`monodromy`) | 4.5 | 16 / 16 | 5 | 7 |

Each scenario comes at two lattice sizes and with or without tensor modes:

| file | lattice | perturbations |
|---|---|---|
| `params_<scenario>_N64.txt` | 64³, small box | scalar |
| `params_<scenario>_N64_tensor.txt` | 64³, small box | scalar + tensor |
| `params_<scenario>_N256.txt` | 256³, paper box | scalar |
| `params_<scenario>_N256_tensor.txt` | 256³, paper box | scalar + tensor |

## Running

Build the example as usual, then run from `Examples/InflatonField`, because
the params files include their spectra table with a path relative to the run
directory (`FILE = paper_examples/spectra/<scenario>_N<N>.txt`):

```bash
cd GRTeclyn/Examples/InflatonField
make -j 8
mpiexec -n 8 ./InflatonField3d.gnu.MPI.ex paper_examples/params_quadratic_N64.txt
```

To run somewhere else, copy or symlink `paper_examples` into the run
directory, or edit the `FILE =` line.

Cost of a full run, measured on Swirles (COSMOS, Cambridge):

| | quadratic | USR | monodromy |
|---|---|---|---|
| steps, 64³ / 256³ | 612 / 1280 | 1019 / 3200 | 1120 / 4480 |
| 64³, 8 MPI ranks × 2 threads (AMD EPYC Genoa) | 3.5 min | 5.5 min | 6 min |
| 64³, 1 Intel PVC GPU (0.1 s per step) | 1.1 min | ≈ 2 min | ≈ 2 min |
| 256³, 1 Intel PVC GPU (2.8 s per step, 22 GB) | ≈ 1 h | ≈ 2.5 h | ≈ 3.5 h |

## What to expect

On these set-ups, the mean field and extrinsic curvature follow the
homogeneous STOIIC background to better than 10⁻⁹ in relative terms, except
in the monodromy case once the field fragments (from about 3 e-folds), and
the runs end at ⟨N⟩ = 7.66, 9.59 and 6.98 e-folds, as in the paper. The
relative Hamiltonian constraint violation (`Ham_absrel`) stays below 10⁻⁶.
The relative momentum constraint violation (`Mom_absrel`) starts at 10⁻⁴ to
10⁻³ with scalar perturbations only. Tensor modes are transverse only up to
the finite-difference error at the cutoff, so the 64³ `_tensor` runs, with
four to six points per cutoff wavelength, start at 15–22% (10⁻³ to 10⁻² at
256³); this decays as the modes leave the Hubble radius. In the monodromy
runs it grows to about 50% once domain walls form, as discussed in the paper.

## Set-up

* Units: the reduced Planck mass is `Mp = 1/sqrt(8 pi G) = 100`
  (`scalar_field.G_Newton = 1/(8 pi 10^4)`).
* Evolution: BSSN without constraint damping, no Kreiss-Oliger dissipation,
  geodesic (synchronous) gauge, fourth-order stencils, `dt = 0.025 dx`,
  periodic box, no AMR.
* Background: STOIIC evolves the homogeneous background and the linear modes
  for ΔN e-folds before the start of the simulation; `init.background_phi`
  and `init.background_dphi` are its values at that time, and GRTeclyn
  recomputes H from the Friedmann equation.
* Scalar perturbations: the comoving curvature mode functions R_k and dR_k from
  STOIIC (`init.use_stoiic_spectra = 1`), drawn with Rayleigh amplitudes and
  uniform phases (`init.random_seed`).
* Tensor perturbations (`*_tensor.txt`): GRTeclyn's analytic Bunch-Davies
  mode functions for the same background. The scalar realisation is the same
  with and without tensors.
* UV cutoff: all modes with k ≥ k_c are removed, with
  k_c = π N_w / (2L) ≈ σ aH. In the small boxes, N_w = N, i.e. the box is the
  largest that resolves the cutoff with four points per wavelength (as in
  ISTORIZ). The paper boxes resolve more than that, so `N_coarse` sets
  N_w = 2σLH/π (122, 81 and 46). The monodromy box is the paper one at both
  sizes.
* Output: about 20 plotfiles per run with all evolved variables, the
  constraints and the `InflationFields` (R, h₊, hₓ). At 256³ each plotfile
  is about 4.5 GB, i.e. about 100 GB per run: restrict `amr.plot_vars` (e.g.
  `amr.plot_vars = chi K phi Pi`) or raise `amr.plot_int` to save disk.
  Checkpoints are off (`amr.check_int = -1`); see the comment in the files
  to turn them on for the long 256³ runs.

## Regenerating the inputs

The files were generated with [STOIIC_GR](https://github.com/YL-codehub/STOIIC_GR)
using the configs in `stoiic/`, through its ISTORIZ writer
(`input_file: "params_file2"`, `stoiic: -1`), and translated to GRTeclyn
parameters by `make_params.py`. For example, for the 64³ quadratic case:

```bash
# in STOIIC_GR, with the configs copied to input/grteclyn/
mpiexec -n 8 python3 main.py --n 64 --nsim 64 --sn 32 --cfg grteclyn/config_quadratic_N64
# back in Examples/InflatonField/paper_examples
python3 make_params.py <STOIIC_GR>/output/dparams.txt quadratic
```

`make_params.py` checks that the Friedmann H of the GRTeclyn potential matches
STOIIC's and that the spectra cover the cutoff. To try other lattices or
boxes, change `scales.length` and `--n/--nsim` (and the cutoff `sigma_R`,
`write.sim.sigma`) in a config.

#!/usr/bin/env python3
"""Translate a STOIIC_GR initial-data file into GRTeclyn InflatonField inputs.

STOIIC_GR (write.sim.software: "ISTORIZ", input_file: "params_file2",
stoiic: -1) writes an ISTORIZ-style dparams.txt: flat GRChombo/ISTORIZ keys
plus the curvature mode functions R_k, dR_k at the start of the simulation.
This script writes, for one scenario and lattice size N,

  spectra/<scenario>_N<N>.txt          the init_k/re_R_k/im_R_k/re_dR_k/im_dR_k
                                       table, read with init.use_stoiic_spectra
  params_<scenario>_N<N>.txt           scalar perturbations only
  params_<scenario>_N<N>_tensor.txt    scalar and tensor perturbations

The params files include the spectra file with the ParmParse "FILE ="
directive, whose path is relative to the run directory (Examples/InflatonField
by default, see --spectra-dir).

Usage:
  python3 make_params.py <dparams.txt> <quadratic|usr|monodromy> [--outdir DIR]
"""

import argparse
import math
import os
import sys

import numpy as np

SPECTRA_KEYS = ["init_k", "re_R_k", "im_R_k", "re_dR_k", "im_dR_k"]

# Scenarios of arXiv:2502.06783. draw_efolds is STOIIC's draw.time, i.e. the
# number of e-folds of linear evolution before the start of the simulation.
SCENARIOS = {
    "quadratic": {
        "title": "quadratic (slow-roll) inflation",
        "potential": "V = m^2 phi^2 / 2",
        "draw_efolds": 5,
        "potential_type": 1,
    },
    "usr": {
        "title": "inflection-point (ultra-slow-roll) inflation",
        "potential": "V = Lambda v^4 (phi^2/3) (3 phi^2 + 2 sqrt(2) phi v + 6 v^2)"
        " / (3 phi^2 + 2 v^2)^2, Lambda = 1.86e-6, v = 0.198 Mp",
        "draw_efolds": 3,
        "potential_type": 8,
    },
    "monodromy": {
        "title": "strong resonance (monodromy) inflation",
        "potential": "V = m^2 phi^2 / 2 + Lambda^4 W(phi - phi_c) [cos((phi - phi_c)/f) - 1],"
        " alpha = 10, b = 1.32, beta = 2",
        "draw_efolds": 5,
        "potential_type": 9,
    },
}


def read_dparams(path):
    """Return (scalars, arrays): scalar keys as strings, spectra as token lists."""
    scalars, arrays = {}, {}
    with open(path) as f:
        for line in f:
            if "=" not in line or line.lstrip().startswith("#"):
                continue
            key, value = line.split("=", 1)
            key = key.strip()
            if key in SPECTRA_KEYS:
                arrays[key] = value.split()
            else:
                scalars[key] = value.split("#", 1)[0].strip()
    missing = [k for k in SPECTRA_KEYS if k not in arrays]
    if missing:
        sys.exit(f"{path}: no {missing}; was it written with input_file: params_file2?")
    lengths = {len(v) for v in arrays.values()}
    if len(lengths) != 1:
        sys.exit(f"{path}: spectra arrays have different lengths {lengths}")
    return scalars, arrays


def potential_keys(ptype, p):
    """GRTeclyn potential keys for a STOIIC/ISTORIZ potential_type and params."""
    if ptype == 1:
        return [("potential.type", "quadratic"), ("scalar_field.scalar_mass", p[1])]
    if ptype == 8:
        return [
            ("potential.type", "usr"),
            ("potential.usr_lambda", p[1]),
            ("potential.usr_v", p[2]),
            # Not used by the usr model, but read by the Potential class
            ("scalar_field.scalar_mass", "0.0"),
        ]
    if ptype == 9:
        return [
            ("potential.type", "monodromy"),
            ("scalar_field.scalar_mass", p[1]),
            ("potential.feature_location", p[2]),
            ("potential.feature_width", p[3]),
            ("potential.feature_amplitude", p[4]),
            ("potential.feature_period", p[5]),
        ]
    sys.exit(f"potential_type {ptype} has no GRTeclyn InflatonField equivalent here")


def potential_value(ptype, p, phi):
    """Python copy of Potential.hpp, used to check the Friedmann background."""
    p = {i: float(v) for i, v in p.items()}
    if ptype == 1:
        return 0.5 * (p[1] * phi) ** 2
    if ptype == 8:
        lam, v = p[1], p[2]
        return (lam * v**4 * phi**2 * (3 * phi**2 + 2 * math.sqrt(2) * phi * v + 6 * v**2)
                / (3 * (3 * phi**2 + 2 * v**2) ** 2))
    if ptype == 9:
        m, loc, width, amp, period = p[1], p[2], p[3], p[4], p[5]
        arg = (phi - loc) / period
        window = 0.25 * (1 + math.tanh(arg)) * (1 + math.tanh((loc - phi + width) / period))
        return 0.5 * (m * phi) ** 2 + amp * window * (math.cos(arg) - 1)
    raise ValueError(ptype)


PARAMS_TEMPLATE = """\
# GRTeclyn InflatonField: {title}
# Example from arXiv:2502.06783, {tensor_desc}.
#   {potential}
#   N = {N}^3 cells, box L = {LH:.6g}/H, UV cutoff k_c = {sigma_eff:.4g} aH
#   (STOIIC sigma = {sigma:g}), {draw_efolds} e-folds of linear evolution
#   before the start, then evolved to t = {duration:.3g}/H (slightly fewer
#   e-folds, as H decreases) in about {n_steps} steps.
# Units: reduced Planck mass Mp = 1/sqrt(8 pi G) = {Mp:.6g}; H = {H:.10g}.
#
# Generated by paper_examples/make_params.py from STOIIC_GR
# (paper_examples/stoiic/config_{scenario}_N{N}.yaml), see paper_examples/README.md.
# Run from Examples/InflatonField so that the spectra file below is found:
#   mpiexec -n 8 ./InflatonField3d.gnu.MPI.ex paper_examples/params_{scenario}_N{N}{suffix}.txt

#################################################
# Filesystem parameters
grteclyn.output_path = ./

amr.verbose = 1

# Ignored by non-GPU builds (ParmParse has no preprocessor, so this line is
# always active)
amrex.use_gpu_aware_mpi = 1

# Plotfiles default to <grteclyn.output_path>/plots and checkpoint files to
# <grteclyn.output_path>/checkpoints
# amr.check_file = checkpoints/chk
# amr.plot_file = plots/plt
# amr.restart = checkpoints/chk00030

# Number of level-0 time steps between checkpoint files; -1 disables output.
# Set e.g. amr.check_int = {plot_int_x5} to be able to restart long runs.
amr.check_int = -1
# Number of level-0 time steps between plot files (about {nplots} per run)
amr.plot_int = {plot_int}

amr.plot_vars = ALL
# R, hplus and hcross (InflationFields) and the Hamiltonian/momentum constraints
amr.derive_plot_vars = constraints InflationFields

#################################################
# Initial Data parameters

# 1 / (8 pi Mp^2), i.e. code units where the reduced Planck mass is {Mp:.6g}
scalar_field.G_Newton = {G_Newton}

# Bunch-Davies perturbations: scalar from the STOIIC spectra, tensor from
# GRTeclyn's analytic mode functions
init.scalar_init = 1
init.tensor_init = {tensor_init}

init.alpha = 0
init.random_seed = {seed}

# Sharp UV cutoff at k_c = pi N_w / (2 L), with N_w = N_coarse if set, else N
init.use_window = 1
init.Delta = 1
{n_coarse_block}
# Homogeneous background at the start of the simulation (from STOIIC)
init.background_phi = {background_phi}
init.background_dphi = {background_dphi}

# Potential: one of quadratic, quadratic_bump, usr, monodromy, punctuated,
# inverted_quadratic_bump, quadratic_step (see Potential.hpp)
{potential_block}

# Scalar mode functions from STOIIC_GR (init_k, re/im_R_k, re/im_dR_k)
init.use_stoiic_spectra = 1
FILE = {spectra_path}

#################################################
# Grid parameters

amr.n_cell = {N} {N} {N}
geometry.prob_extent = {L} {L} {L}

amr.max_level = 0 # There are (max_level+1) grids, so min is zero

# Regridding interval in level timesteps; -1 disables regridding
amr.regrid_int = -1

# Maximum number of cells along each side of a grid box
amr.max_grid_size = {max_grid_size}
# Grid box side lengths must be divisible by this value
amr.blocking_factor = {blocking_factor}

# Number of cells by which tagged regions are grown
amr.n_error_buf = 0

#################################################
# Boundary condition parameters

# Periodic directions: 0 = false, 1 = true
geometry.is_periodic = 1 1 1

#################################################
# Evolution parameters

# dt is dx * dt_multiplier on each grid level
evolution.dt_multiplier = {dt_multiplier}
# {duration:.3g}/H of cosmic time
evolution.stop_time = {stop_time}
# evolution.max_steps = 1000000

# Spatial derivative order (only affects the CCZ4 RHS)
evolution.spatial_derivative_order = 4

evolution.nan_check = true

# Geodesic (synchronous) gauge: lapse = 1, shift = 0
gauge.lapse_advec_coeff = 0.0
gauge.lapse_coeff = 0.0
gauge.lapse_power = 1.0
gauge.shift_advec_coeff = 0.0
gauge.shift_Gamma_coeff = 0.0
gauge.eta = 0.0

# BSSN without constraint damping
ccz4.formulation = 1 # 1 for BSSN, 0 for CCZ4
ccz4.kappa1 = 0.
ccz4.kappa2 = 0.
ccz4.kappa3 = 0.
ccz4.covariantZ4 = 0 # 0: keep kappa1; 1 [default]: replace kappa1 -> kappa1/lapse

# coefficient for KO numerical dissipation
evolution.sigma = 0.0

ccz4.min_chi = 1.e-30
ccz4.min_lapse = 1.e-30
"""


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("dparams", help="dparams.txt written by STOIIC_GR")
    parser.add_argument("scenario", choices=SCENARIOS)
    parser.add_argument("--outdir", default=os.path.dirname(os.path.abspath(__file__)),
                        help="where to write the files (default: next to this script)")
    parser.add_argument("--spectra-dir", default="paper_examples/spectra",
                        help="spectra directory as seen from the run directory")
    parser.add_argument("--nplots", type=int, default=20,
                        help="approximate number of plotfiles per run")
    parser.add_argument("--seed", type=int, default=3539263, help="init.random_seed")
    args = parser.parse_args()

    scen = SCENARIOS[args.scenario]
    s, arrays = read_dparams(args.dparams)

    ptype = int(s["potential_type"])
    if ptype != scen["potential_type"]:
        sys.exit(f"potential_type {ptype} does not match scenario {args.scenario}")
    p = {i: s[f"potential_param_{i}"] for i in range(1, 6)}

    G_Newton = float(s["G_Newton"])
    Mp = 1 / math.sqrt(8 * math.pi * G_Newton)
    phi0, dphi0 = float(s["background_phi"]), float(s["background_dphi"])
    H = float(s["hubble"])
    N, L = int(s["N"]), float(s["L"])
    sigma = float(s["sto_sigma"])
    dt_multiplier, stop_time = float(s["dt_multiplier"]), float(s["stop_time"])

    # The background GRTeclyn computes from Friedmann must be STOIIC's
    V0 = potential_value(ptype, p, phi0)
    H_friedmann = math.sqrt(8 * math.pi * G_Newton / 3 * (0.5 * dphi0**2 + V0))
    if abs(H_friedmann / H - 1) > 1e-6:
        sys.exit(f"Friedmann H = {H_friedmann} differs from STOIIC hubble = {H}")

    # GRTeclyn's window cuts at pi N_w / (2L): N_w = N_coarse puts it at sigma aH
    # when the grid can resolve more than that.
    n_coarse = round(2 * sigma * L * H / math.pi)
    if n_coarse < N:
        n_w = n_coarse
        n_coarse_block = (
            "# Cutoff at sigma aH = pi N_coarse / (2L) rather than at the grid scale\n"
            "# (read from the top-level table)\n"
            f"N_coarse = {n_coarse}\n")
    else:
        n_w = N
        n_coarse_block = ""
    k_c = math.pi * n_w / (2 * L)
    k_table = np.array(arrays["init_k"], dtype=float)
    if k_table.max() < k_c:
        sys.exit(f"spectra stop at k = {k_table.max()} < cutoff {k_c}: increase sigma_R")

    n_steps = round(stop_time / (dt_multiplier * L / N))
    plot_int = max(1, round(n_steps / args.nplots))
    max_grid_size, blocking_factor = (32, 16) if N <= 64 else (N // 2, 32)

    spectra_name = f"{args.scenario}_N{N}.txt"
    os.makedirs(os.path.join(args.outdir, "spectra"), exist_ok=True)
    with open(os.path.join(args.outdir, "spectra", spectra_name), "w") as f:
        f.write(f"# STOIIC_GR mode functions for params_{args.scenario}_N{N}*.txt:\n"
                "# comoving curvature perturbation R_k and its cosmic-time derivative\n"
                f"# dR_k at the start of the simulation, on {len(k_table)} k values.\n")
        for key in SPECTRA_KEYS:
            f.write(f"{key} = " + " ".join(arrays[key]) + "\n")

    fields = dict(
        title=scen["title"], potential=scen["potential"], scenario=args.scenario,
        N=N, L=s["L"], LH=L * H, sigma=sigma, sigma_eff=k_c / H,
        draw_efolds=scen["draw_efolds"], duration=stop_time * H, n_steps=n_steps,
        Mp=Mp, H=H, plot_int=plot_int, plot_int_x5=5 * plot_int, nplots=args.nplots,
        G_Newton=s["G_Newton"], seed=args.seed, n_coarse_block=n_coarse_block,
        background_phi=s["background_phi"], background_dphi=s["background_dphi"],
        potential_block="\n".join(f"{k} = {v}" for k, v in potential_keys(ptype, p)),
        spectra_path=f"{args.spectra_dir}/{spectra_name}",
        max_grid_size=max_grid_size, blocking_factor=blocking_factor,
        dt_multiplier=s["dt_multiplier"], stop_time=s["stop_time"],
    )
    for tensor_init, suffix, desc in [(0, "", "scalar perturbations only"),
                                      (1, "_tensor", "scalar and tensor perturbations")]:
        name = f"params_{args.scenario}_N{N}{suffix}.txt"
        with open(os.path.join(args.outdir, name), "w") as f:
            f.write(PARAMS_TEMPLATE.format(tensor_init=tensor_init, suffix=suffix,
                                           tensor_desc=desc, **fields))
        print(f"wrote {name}")
    print(f"wrote spectra/{spectra_name}: {len(k_table)} k values, cutoff "
          f"k_c = {k_c / H:.4g} aH (N_w = {n_w}), L = {L * H:.6g}/H, "
          f"{n_steps} steps, plot_int = {plot_int}, H/H_friedmann - 1 = {H / H_friedmann - 1:.2e}")


if __name__ == "__main__":
    main()

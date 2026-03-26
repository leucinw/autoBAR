# autoBAR: Free Energy Simulator with Tinker Software Packages

[![Python 3.8+](https://img.shields.io/badge/python-3.8%2B-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)
[![Free Energy](https://img.shields.io/badge/method-BAR-orange.svg)](https://en.wikipedia.org/wiki/Bennett_acceptance_ratio)

A lightweight automation tool for alchemical free energy simulations using the
[Bennett Acceptance Ratio (BAR)](https://en.wikipedia.org/wiki/Bennett_acceptance_ratio)
method with the polarizable **AMOEBA** and **AMOEBA+** force fields.

BAR is a statistically optimal estimator of the free energy difference between
two states — it uses overlap information from both endpoints to minimize variance.

---

## How it works

```
setup  →  dynamic  →  bar  →  result
  │           │          │        │
Generate   Run MD    Run BAR   Parse .ene
key/xyz   at each   analysis   & write
files     λ window  on .arc    result.txt
```

---

## Prerequisites

### Python packages

| Package | Version |
|---------|---------|
| `ruamel.yaml` | `<0.18` |
| `numpy` | any recent |
| `scipy` | `>=1.16` |

### Compiled binaries

| Binary | Notes |
|--------|-------|
| **Tinker** (CPU) | Required for gas-phase simulations |
| **Tinker-GPU** | Required for liquid-phase simulations |

---

## Installation / Quick-start

```shell
# 1. Clone the repository
git clone https://github.com/leucinw/autoBAR.git
cd autoBAR

# 2. Install Python dependencies
pip install "ruamel.yaml<0.18" numpy "scipy>=1.16"
```

---

## Required input files

| File | Description |
|------|-------------|
| `gas_xyz` | Ligand/small-molecule Tinker XYZ file |
| `box_xyz` | Simulation box Tinker XYZ file (with box info in the second line) |
| `parameters` | Tinker force-field parameter file |
| `settings.yaml` | autoBAR run settings (see [Configuration reference](#configuration-reference)) |

> **Note:** User-customized Tinker key files are supported — see [Notes / Advanced](#notes--advanced).

---

## Usage

### Interactive mode

Run each step individually as needed:

```shell
# Step 1 — generate Tinker input files (XYZ symlinks, KEY files)
python autoBAR.py setup

# Step 2 — submit MD simulations at every λ window
python autoBAR.py dynamic

# Step 3 — run BAR analysis once MD trajectories are complete
python autoBAR.py bar

# Step 4 — summarize results; writes result.txt
python autoBAR.py result
```

### Automated mode

```shell
# Runs all four steps sequentially, polling until each step is complete
python autoBAR.py auto
```

Each step produces:

| Step | Output |
|------|--------|
| `setup` | `{phase}/` directory with `.xyz`, `.key` files |
| `dynamic` | `{phase}/*.arc` trajectory files |
| `bar` | `{phase}/*.ene` free-energy files |
| `result` | `result.txt` with free energies per window and total ΔG |

---

## Examples

Two example systems are provided in the `examples/` directory.

### Ion HFE (liquid phase only)

```
├── liquid/
├── Na-water.xyz
├── Na.xyz
├── result.txt          ← final free energy (ΔG in kcal/mol)
├── settings.yaml
└── water03.prm
```

### Phenol HFE (gas + liquid phases)

```
├── amoeba09.prm
├── gas/
├── liquid/
├── phenol_solv.xyz
├── phenol.xyz
├── result.txt          ← gas + liquid contributions and total ΔG
└── settings.yaml
```

`result.txt` contains the free energy contribution of every λ → λ′ window,
followed by the total ΔG (sum over all windows) and, if one-step FEP was used,
the perturbed free energies.

---

## Configuration reference

The following table covers the most commonly used `settings.yaml` keys.
See [`dat/settings.yaml`](dat/settings.yaml) for the full annotated example.

| Key | Type | Default | Description |
|-----|------|---------|-------------|
| `gas_xyz` | string | — | Ligand Tinker XYZ file (required) |
| `box_xyz` | string | — | Simulation box Tinker XYZ file (required) |
| `parameters` | string | — | Tinker parameter file (required) |
| `lambda_window` | string | `"default"` | Lambda spacing: `"coarser"` (fewer windows) or `"default"` (standard) |
| `checking_time` | float | `120.0` | Polling interval (s) used in `auto` mode |
| `verbose` | int | `1` | Set to `0` to suppress progress output |
| `copy_arc_for_perturb` | bool | `false` | Reuse end-state `.arc` for FEP windows instead of re-running MD |
| `liquid_md_total_time` | float | — | Liquid MD simulation length (ns) |
| `liquid_md_time_step` | float | — | Integration timestep for liquid MD (fs) |
| `liquid_md_write_freq` | float | — | Snapshot write frequency for liquid MD (ps) |
| `liquid_md_ensemble` | string | `"NPT"` | Liquid MD ensemble (`"NPT"` or `"NVT"`) |
| `liquid_md_temperature` | float | — | Liquid MD temperature (K) |
| `liquid_md_pressure` | float | — | Liquid MD pressure (atm) |
| `gas_md_total_time` | float | — | Gas MD simulation length (ns); set to `0` to skip gas phase |
| `gas_md_time_step` | float | — | Integration timestep for gas MD (fs) |
| `gas_md_write_freq` | float | — | Snapshot write frequency for gas MD (ps) |
| `gas_md_temperature` | float | — | Gas MD temperature (K) |
| `node_list` | list | — | Compute node names (required unless running as developer) |
| `liquid_key` | string | built-in | Custom Tinker key file for the liquid phase |
| `gas_key` | string | built-in | Custom Tinker key file for the gas phase |
| `lambda_window_file` | string | — | Path to a custom λ-window file (overrides `lambda_window`) |
| `manual_ele_scale` | bool | `false` | Enable manual electrostatic scaling (see Notes) |

---

## Notes / Advanced

### Minimal settings for HFE simulations

```yaml
lambda_window: coarser      # Fewer windows; minimal accuracy loss
liquid_md_total_time: 1.25  # 1.25 ns total (last 4/5 = 1 ns used in BAR)
liquid_md_time_step: 2.0    # fs; works well with RESPA integrator
gas_md_total_time: 1.25     # 1.25 ns total
gas_md_time_step: 0.1       # fs; small step for stochastic gas-phase dynamics
```

### One-step perturbation (FEP)

Place a file named `{parameters}.prm_XX` (e.g., `final.prm_01`) with a small
parameter perturbation in the working directory.  autoBAR will automatically
detect it and treat it as an extra end-state.

- Multiple `*.prm_XX` files are supported (up to 99 simultaneous perturbations).
- No changes to `settings.yaml` are needed.

### Manual electrostatic scaling

Set `manual_ele_scale: true` in `settings.yaml`.  The electrostatic parameters
will be scaled explicitly in each `.key` file using `utils/elescale.py`.

### Custom key files

Provide `liquid_key` and/or `gas_key` in `settings.yaml` to override the
built-in Tinker key templates.  You are responsible for the correctness of
custom key files.

### Lambda-window spelling note

Both `coarser` and `courser` (legacy spelling) are accepted for `lambda_window`.
The built-in coarser window set is stored as `dat/orderparams_courser` and
`dat/orderparams_coarser` (identical content).

---

## Citation / Contact

If you use autoBAR in your research, please cite the relevant Tinker and AMOEBA
publications and acknowledge this tool.

**Author:** Chengwen Liu  
**Email:** liuchw2010@gmail.com  
**Affiliation:** University of Texas at Austin

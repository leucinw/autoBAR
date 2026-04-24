
# autoBAR

**Automated Free Energy Simulations with the Bennett Acceptance Ratio (BAR) Method**

autoBAR is a lightweight automation tool for running alchemical free energy simulations using [Tinker](https://dasher.wustl.edu/tinker/) and [Tinker9](https://github.com/TinkerTools/tinker9). It supports the polarizable AMOEBA and AMOEBA+ force fields and uses the [Bennett Acceptance Ratio](https://en.wikipedia.org/wiki/Bennett_acceptance_ratio) for free energy estimation.

---

## Prerequisites

| Category | Requirements |
|----------|-------------|
| **Python** (≥ 3.8) | `numpy`, `scipy ≥ 1.16`, `ruamel.yaml < 0.18` |
| **Simulation Software** | [Tinker](https://dasher.wustl.edu/tinker/) (CPU) and [Tinker9](https://github.com/TinkerTools/tinker9) (GPU) |

## Project Structure

```
autoBAR/
├── autoBAR.py              # Main driver script
├── dat/                    # Default configuration files
│   ├── settings.yaml       # Example settings template
│   ├── gas.key             # Default gas-phase key file
│   ├── liquid.key          # Default liquid-phase key file
│   ├── orderparams_courser # Coarser lambda schedule (17 windows)
│   ├── orderparams_default # Standard lambda schedule (26 windows)
│   └── tinker.env          # Tinker executable paths
├── utils/                  # Utility modules
│   ├── checkautobar.py     # Progress monitoring for MD and BAR jobs
│   ├── elescale.py         # Electrostatic parameter scaling
│   ├── parmOPT.py          # Parameter optimization targeting experimental HFE and liquid density
│   └── submitTinker.py     # Cluster job submission helper
└── examples/               # Ready-to-use example systems
    ├── Ion-HFE/            # Na⁺ hydration free energy
    └── Phenol-HFE/         # Phenol hydration free energy
```

## Setup

### 1. Configure Tinker Paths

Edit `dat/tinker.env` to point to your local Tinker and Tinker9 installations:

```bash
export TINKER8=/path/to/tinker
export DYNAMIC8="$TINKER8/dynamic"
export BAR8="$TINKER8/bar"
export MINIMIZE8="$TINKER8/minimize"

export tk9home=/path/to/tinker9/build
export DYNAMIC9="$tk9home/dynamic9"
export BAR9="$tk9home/bar9"
export MINIMIZE9="$tk9home/minimize9"
```

### 2. Prepare Input Files

Place the following four files in your working directory:

| File | Description |
|------|-------------|
| `gas_xyz` | Ligand Tinker `.xyz` file |
| `box_xyz` | Solvated system Tinker `.xyz` file (box dimensions on line 2) |
| `parameters` | Tinker force field parameter file (`.prm`) |
| `settings.yaml` | Simulation settings — see the [template](dat/settings.yaml) |

> **Tip:** Custom `.key` files are supported. Set `liquid_key` and/or `gas_key` in `settings.yaml` to use your own key files instead of the defaults.

### 3. Configure `settings.yaml`

Key settings to configure (see [`dat/settings.yaml`](dat/settings.yaml) for the full reference):

```yaml
# Lambda schedule: "courser" (17 windows) or "default" (26 windows)
lambda_window: courser

# Input files
gas_xyz: ligand.xyz
box_xyz: solvated.xyz
parameters: forcefield.prm

# Liquid-phase MD settings
liquid_md_total_time: 1.25    # ns
liquid_md_time_step: 2.0      # fs (RESPA integrator)
liquid_md_ensemble: NPT
liquid_md_temperature: 300.0  # K

# Gas-phase MD settings (set gas_md_total_time to 0 to skip)
gas_md_total_time: 1.25       # ns
gas_md_time_step: 0.1         # fs (stochastic dynamics)
gas_md_temperature: 300.0     # K

# Cluster nodes for job distribution
node_list:
  - node01
  - node02
```

## Usage

autoBAR can be run in **interactive** or **automated** mode.

### Interactive Mode

Run individual steps as needed:

```bash
# Step 1: Generate simulation input files
python autoBAR.py setup

# Step 2: Run MD simulations across all lambda windows
python autoBAR.py dynamic

# Step 3: Perform BAR free energy analysis
python autoBAR.py bar

# Step 4: Collect and summarize results
python autoBAR.py result
```

### Automated Mode

Run the entire workflow end-to-end:

```bash
python autoBAR.py auto
```

In automated mode, autoBAR will:
1. Generate all input files (`setup`)
2. Submit MD jobs and wait for completion (`dynamic`)
3. Submit BAR analysis jobs and wait for completion (`bar`)
4. Collect and print the final free energy (`result`)

## Workflow Overview

```
setup ──► dynamic ──► bar ──► result
  │          │         │        │
  │          │         │        └─ Parse .ene files, sum ΔG, write result.txt
  │          │         └────────── Run BAR on trajectory pairs (arc0, arc1)
  │          └──────────────────── Run MD at each λ window (liquid: GPU, gas: CPU)
  └─────────────────────────────── Generate .key/.xyz per λ window, minimize
```

## Examples

### Ion Hydration Free Energy (Na⁺)

Only the liquid phase is needed (small ion → gas phase is skipped automatically):

```
examples/Ion-HFE/
├── Na.xyz              # Gas-phase ion structure
├── Na-water.xyz        # Solvated ion + water box
├── water03.prm         # Force field parameters
├── settings.yaml       # Simulation settings
└── result.txt          # Reference output
```

### Phenol Hydration Free Energy

Both gas and liquid phases are required:

```
examples/Phenol-HFE/
├── phenol.xyz          # Gas-phase phenol structure
├── phenol_solv.xyz     # Solvated phenol + water box
├── amoeba09.prm        # AMOEBA force field parameters
├── settings.yaml       # Simulation settings
└── result.txt          # Reference output
```

## One-Step Free Energy Perturbation (FEP)

autoBAR supports one-step perturbation for evaluating small parameter changes:

1. Place one or more `{parameters}.prm_XX` files (where `XX` = `01` to `99`) in the working directory alongside your main `.prm` file
2. No changes to `settings.yaml` are needed — autoBAR detects the files automatically
3. Each perturbed parameter file defines a new end state
4. Results are reported as `FEP_001`, `FEP_002`, etc. in `result.txt`

This is useful for sensitivity analysis or parameter optimization (see `utils/parmOPT.py`).

## Parameter Optimization

The `utils/parmOPT.py` utility optimizes force field parameters to simultaneously match experimental hydration free energy (HFE) and neat liquid density at one or more temperatures:

```bash
python utils/parmOPT.py
```

It uses `scipy.optimize.least_squares` (TRF, soft-L1 loss) with a custom Jacobian. The HFE row is evaluated by one-step FEP via autoBAR; each density row uses the fluctuation formula (Wang et al. 2013, Eq. 4) applied to the most recent per-temperature trajectory, with all `$ANALYZE` calls dispatched in parallel alongside the autoBAR HFE run.

### Required settings

Add these fields to `settings.yaml` (copy `utils/settings.yaml` as a starting point):

```yaml
# --- parmOPT required ---
parameters: forcefield.prm               # Tinker parameter file
expt_hfe: -4.75                          # Experimental HFE (kcal/mol)

# Single-temperature density target (scalar form):
temperature: 298.15                      # Simulation temperature (K)
expt_density: 997.0                      # Experimental density in kg/m³
                                         # NOTE: units are kg/m³ — water ≈ 997, ethanol ≈ 789

liquid_dir: neat_liquid                  # Directory for neat-liquid MD (must NOT be "liquid"
                                         # — autoBAR already uses ./liquid/ for HFE windows)
n_production: 500                        # Number of production frames per MD run

opt_params: "vdwpair-401-402 3.8 0.05"  # Force field term + initial parameter values
params_range: "0.5 0.02"                # Search range (±) for each parameter

# --- parmOPT optional ---
hfe_weight: 1.0                          # Weight applied to the HFE residual (default: 1.0)
density_weight: 1.0                      # Weight applied to the density residual (default: 1.0)
liquid_base: neat_liq                    # Coordinate/trajectory base name inside liquid_dir
                                         # (default: "neat_liq"); sets xyz/key/sh/arc/dyn names
liquid_key: neat_liq                     # Key file basename (default: same as liquid_base)
n_equil: 200                             # Frames to discard as equilibration (default: 200)
```

### Choosing the weights

The optimizer minimizes the sum of squared weighted residuals:

```
cost = (hfe_weight × ΔHFE)² + Σ (density_weight_i × Δρ_i)²
```

Because HFE is in kcal/mol and density is in kg/m³, a raw error of 1 kcal/mol
and a raw error of 1 kg/m³ are physically very different but contribute equally to
the cost when both weights are 1.0. This imbalance will cause the optimizer to
focus almost entirely on density and neglect HFE. The weights should be chosen so
that a *chemically meaningful* error in each property produces a similar cost
contribution.

**Rule of thumb — inverse acceptable error:**

| Property | Typical acceptable error | Suggested weight |
|----------|--------------------------|------------------|
| HFE | 0.2 kcal/mol | `hfe_weight: 5.0` |
| Density (kg/m³) | 2 kg/m³ | `density_weight: 0.5` |

Both give a weighted residual of ~1 at the acceptable-error threshold, so neither
property dominates the other.

You can also set the weights relative to each other without fixing a scale.
A useful starting point is:

```
density_weight / hfe_weight ≈ (target HFE precision) / (target density precision)
                             = 0.2 kcal/mol / 2 kg/m³  =  0.1
```

For example, `hfe_weight: 1.0` and `density_weight: 0.1` means a 1 kcal/mol HFE
error costs the same as a 10 kg/m³ density error.

**Practical starting suggestion for water-like solvents:**

```yaml
hfe_weight: 1.0
density_weight: 0.1
```

If the optimized parameters drift too far from the experimental density while
fitting HFE well (or vice versa), increase the weight of the lagging property by
a factor of 2–5 and rerun. Monitor both residuals in `parmOPT.log` at each step
to guide the adjustment.

### Multi-temperature density fitting

Providing densities at multiple temperatures simultaneously constrains the parameter search and reduces overfitting.

Replace the scalar `temperature` / `expt_density` keys with their list counterparts:

```yaml
# Replace these scalar keys:
#   temperature: 298.15
#   expt_density: 997.0

# With these list keys (lengths must match):
temperatures:    [278.15, 298.15, 318.15]  # temperatures in K
expt_densities:  [999.9,  997.0,  989.3]   # experimental densities in kg/m³ at each T
density_weights: [1.0,    1.0,    1.0]      # optional per-T weights; replaces density_weight
```

> **Important:** `expt_densities` values must be in **kg/m³**, not g/cm³.
> Common reference values: water 278 K → 999.9, 298 K → 997.0, 318 K → 989.3 kg/m³.

**Weight normalization:** each `density_weight` is divided by the number of temperatures internally, so the total density contribution to the cost is the same whether you fit at one temperature or ten. The user-specified weights control the *relative* importance of each temperature; the *aggregate* HFE/density balance is unchanged. The effective weights are printed to `parmOPT.log` at startup.

autoBAR runs the HFE calculation once per optimizer call. Liquid MD is run sequentially for each temperature, producing temperature-tagged arc files (`neat_liq_298.2K.arc`, etc.) so that dyn-file checkpointing and parallel `$ANALYZE` calls work correctly across all temperatures.

### Neat-liquid directory layout

> **Important:** set `liquid_dir` to a name other than `liquid` (e.g., `neat_liquid`).
> autoBAR creates `./liquid/` and `./gas/` for the HFE alchemical windows; a conflicting
> `liquid_dir: liquid` would mix those files with the pure-liquid MD trajectories.

`liquid_dir` must contain files named after `liquid_base` (default `neat_liq`):

| File | Description |
|------|-------------|
| `neat_liq.xyz` | Simulation box Tinker coordinates |
| `neat_liq.key` (or the name set by `liquid_key`) | Tinker key file |
| `neat_liq.sh` | Shell script with a `$DYNAMIC neat_liq ...` line — parsed for timestep, write interval, integrator type, and pressure |

Example directory after a run:
```
neat_liquid/
├── neat_liq.xyz
├── neat_liq.key
├── neat_liq.sh
├── neat_liq.arc          # trajectory (single temperature)
└── neat_liq.dyn          # Tinker checkpoint for dyn-file reuse
```

For multi-temperature runs the arc and dyn files are tagged by temperature:
```
neat_liquid/
├── neat_liq_278.2K.arc
├── neat_liq_278.2K.dyn
├── neat_liq_298.2K.arc
└── neat_liq_298.2K.dyn
```

## Recommended Minimal Settings

For hydration free energy calculations:

| Setting | Value | Notes |
|---------|-------|-------|
| `lambda_window` | `courser` | 17 windows — good accuracy with less cost |
| `liquid_md_total_time` | 1.25 ns | Last 80% (1 ns) used in BAR analysis |
| `liquid_md_time_step` | 2.0 fs | Works well with RESPA integrator |
| `gas_md_total_time` | 1.25 ns | Last 80% (1 ns) used in BAR analysis |
| `gas_md_time_step` | 0.1 fs | Required for gas-phase stochastic dynamics |

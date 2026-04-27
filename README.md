
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
export ANALYZE8="$TINKER8/analyze"
export BAR8="$TINKER8/bar"
export MINIMIZE8="$TINKER8/minimize"

export tk9home=/path/to/tinker9/build
export DYNAMIC9="$tk9home/dynamic9"
export ANALYZE9="$tk9home/analyze9"
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

It uses `scipy.optimize.least_squares` (TRF, soft-L1 loss) with a custom Jacobian. The HFE row is evaluated by one-step FEP via autoBAR; each density row uses the fluctuation formula (Wang et al. 2013, Eq. 4) applied to the most recent per-temperature trajectory, with all `$ANALYZE9` jobs submitted to the GPU cluster in parallel alongside the autoBAR HFE run.

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
production_time: 2.0                     # Production simulation time (ns)

opt_params: "vdwpair-401-402 3.8 0.05"  # Force field term + initial parameter values
params_range: "0.5 0.02"                # Search range (±) per parameter; use 0 to fix

# --- Liquid MD settings (shared between HFE and neat-liquid simulations) ---
liquid_md_time_step:  2.0               # Integration timestep (fs)
liquid_md_write_freq: 0.1               # Trajectory output interval (ps)
liquid_md_pressure:   1.0               # Pressure (atm)

# --- parmOPT optional ---
hfe_weight: 1.0                          # Relative importance weight for the HFE residual (default: 1.0)
density_weight: 1.0                      # Relative importance weight for the density residual (default: 1.0)
hfe_denom: 2.24                          # Scale denominator for HFE (kcal/mol).
                                         # Default: sqrt(|expt_hfe|). Override to fix the normalization
                                         # scale across multiple optimization runs.
density_denom: 31.6                      # Scale denominator for density (kg/m³).
                                         # Default: std-dev of expt_densities (multi-T) or
                                         # sqrt(|expt_density|) (single-T). Override as needed.
liquid_base: neat_liq                    # Coordinate/trajectory base name inside liquid_dir
                                         # (default: "neat_liq"); sets xyz/key/sh/arc/dyn names
liquid_key: neat_liq                     # Key file basename (default: same as liquid_base)
equil_time: 0.5                          # Equilibration time (ns) — discarded before averaging (default: 0.02)
hfe_liquid_key: path/to/liquid.key       # Override the HFE liquid key template used to auto-generate neat_liq.key
                                         # (default: dat/liquid.key in the autoBAR repo)
checking_time: 60                        # Polling interval (s) for MD and analyze job completion (default: 60)
```

### Weights and scale normalization

The optimizer minimizes:

```
cost = (hfe_weight × ΔHFE / hfe_denom)² + Σ (density_weight_i × Δρ_i / density_denom)²
```

There are two separate knobs:

| Knob | Purpose | Settings key |
|------|---------|--------------|
| **Denom** | Normalizes the physical scale of each property so residuals are dimensionless | `hfe_denom`, `density_denom` |
| **Weight** | Controls the *relative importance* of HFE vs. density after normalization | `hfe_weight`, `density_weight` |

#### Scale denominators (`hfe_denom`, `density_denom`)

The denominator converts each raw residual (in its physical unit) into a dimensionless number. parmOPT computes sensible defaults automatically — the same convention used by ForceBalance:

| Property | Default `Denom` | Rationale |
|----------|-----------------|-----------|
| HFE | `sqrt(\|expt_hfe\|)` | e.g. expt\_hfe = −5 kcal/mol → denom ≈ 2.24 kcal/mol |
| Density (single T) | `sqrt(\|expt_density\|)` | e.g. expt\_density = 997 kg/m³ → denom ≈ 31.6 kg/m³ |
| Density (multi T) | `std_dev(expt_densities)` | spread of the experimental values across temperatures |

With these defaults a typical HFE error of ~1 kcal/mol and a typical density error of ~10 kg/m³ both produce a normalized residual of ~0.3–0.5, so `hfe_weight = density_weight = 1.0` gives a balanced starting point **without any manual tuning**.

Override the defaults only when you want to lock the normalization scale across multiple runs (e.g. to make cost values comparable between different optimization attempts):

```yaml
hfe_denom: 2.24      # kcal/mol — fix at sqrt(5) regardless of expt_hfe
density_denom: 31.6  # kg/m³   — fix at sqrt(997) regardless of expt_density
```

#### Relative importance weights (`hfe_weight`, `density_weight`)

Once the denominators have put both residuals on the same scale, the weights express how much you care about one property relative to the other. With the default denominators, `hfe_weight = density_weight = 1.0` is a reasonable starting point for most organic solvents.

Adjust when the optimization consistently sacrifices one property for the other:

- If the optimizer fits HFE well but density drifts — increase `density_weight` (or decrease `hfe_weight`).
- If density is well reproduced but HFE error is large — increase `hfe_weight`.
- A factor-of-2 change in a weight shifts the cost ratio by 4×; start with small adjustments.

The `WtNormRes` column in `parmOPT.log` shows `weight × Δ / denom` for each property at every optimizer step. Watch this column to see which property is driving the cost.

**Practical starting point (works for most solvents):**

```yaml
# Let denom defaults handle unit normalization; use equal weights.
hfe_weight: 1.0
density_weight: 1.0
```

If you have a strong prior on acceptable errors you can also set weights as inverse acceptable errors *after* normalization. For example, if you want a 0.5 kcal/mol HFE error to cost the same as a 5 kg/m³ density error, and the defaults give `hfe_denom ≈ 2.24` and `density_denom ≈ 31.6`:

```
# normalized acceptable errors
hfe_threshold   = 0.5 / 2.24 ≈ 0.22
density_threshold = 5 / 31.6 ≈ 0.16

# set weights so both thresholds give cost contribution ≈ 1
hfe_weight     = 1 / 0.22 ≈ 4.5   → round to 5.0
density_weight = 1 / 0.16 ≈ 6.3   → round to 6.0
```

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

autoBAR runs the HFE calculation once per optimizer call. Neat-liquid MD jobs are submitted to the GPU cluster in parallel with the HFE run: each temperature gets its own coordinate symlink, run script, and GPU card, so all temperatures run simultaneously. Before submitting, parmOPT checks the existing `.arc` file for each temperature — if it already contains the required number of frames (e.g., after a restart), the submission is skipped; if the run is partial and a `.dyn` checkpoint exists, a resume script is generated with only the remaining steps. parmOPT waits for the HFE job and all neat-liquid MD jobs to finish before moving to the next optimizer step.

For the Jacobian, parmOPT first strips the equilibration frames from each trajectory, writing a `*-prod.arc` file containing only production snapshots. It then writes one `$ANALYZE9` run script per (parameter perturbation, temperature) pair, submits all scripts to the GPU cluster simultaneously via `submitTinker.py`, and waits for all output logs before assembling the density sensitivity rows of the Jacobian.

### Optimizing multiple parameter groups

`opt_params` and `params_range` accept either a single string (one group) or a YAML list (multiple groups). Each list entry is one independent force-field term; all groups are optimized simultaneously.

```yaml
opt_params:
  - "vdw-401    3.80 0.050"   # R and epsilon on atom type 401
  - "vdw-402    3.60 0.060"   # R and epsilon on atom type 402
  - "chgpen-403 5.00 0.800"   # two chgpen parameters on atom type 403
params_range:
  - "0.30 0.02"
  - "0.30 0.02"
  - "0.50 0.05"
```

The optimizer sees a single flat parameter vector formed by concatenating the free parameters from all groups; `write_prm` splits it back and writes one line per group.

### Fixing individual parameters

Use `0` in `params_range` to hold a parameter constant. That parameter is excluded from the optimizer but still written to the parameter file at its initial value:

```yaml
# Optimize R only; keep epsilon fixed at 0.050
opt_params:   "vdw-401 3.80 0.050"
params_range: "0.30 0"
```

This works for both the single-string and list forms:

```yaml
opt_params:
  - "vdw-401    3.80 0.050"   # optimize R only
  - "chgpen-403 5.00 0.800"   # optimize both params
params_range:
  - "0.30 0"      # 0 → epsilon is fixed
  - "0.50 0.05"
```

Fixed parameters are annotated as `(fixed)` in each step's log line in `parmOPT.log`.

### Neat-liquid directory layout

> **Important:** set `liquid_dir` to a name other than `liquid` (e.g., `neat_liquid`).
> autoBAR creates `./liquid/` and `./gas/` for the HFE alchemical windows; a conflicting
> `liquid_dir: liquid` would mix those files with the pure-liquid MD trajectories.

The only file you need to supply in `liquid_dir` is the Tinker coordinate file.
The key and per-temperature shell scripts are auto-generated at startup:

| File | Source | Description |
|------|--------|-------------|
| `neat_liq.xyz` | **User-supplied** | Simulation box Tinker coordinates |
| `neat_liq.key` | Auto-generated (if absent) | Derived from the HFE liquid key template by removing `vdw-annihilate`, `vdw-lambda`, `ele-lambda`, and `ligand` lines; shared by all temperatures |
| `neat_liq_{T}K.xyz` | Auto-generated (symlink) | Per-temperature symlink → `neat_liq.xyz`; Tinker writes `neat_liq_{T}K.arc` so parallel GPU runs don't overwrite each other |
| `neat_liq_{T}K.key` | Auto-generated (symlink) | Per-temperature symlink → `neat_liq.key`; ensures Tinker names its output after the temperature-tagged coordinate file |
| `neat_liq_{T}K.sh` | Auto-generated (always) | Per-temperature NPT run script (sources `dat/tinker.env`, calls `$DYNAMIC9`); each is submitted to a separate GPU card |

Example directory after a two-temperature run at 278 K and 298 K with Jacobian computed:
```
neat_liquid/
├── neat_liq.xyz                        # user-supplied
├── neat_liq.key                        # shared key file (PARAMETERS updated each opt step)
├── neat_liq_278K.xyz                   # symlink → neat_liq.xyz
├── neat_liq_278K.key                   # symlink → neat_liq.key
├── neat_liq_278K.sh                    # MD run script for 278 K
├── neat_liq_278K-md.log                # MD log
├── neat_liq_278K.arc                   # full trajectory (equil + production)
├── neat_liq_278K.dyn                   # MD checkpoint for resume
├── neat_liq_278K-prod.arc              # production-only arc (equil frames stripped)
├── neat_liq_278K-prm01.key             # analyze key: prm perturb +Δ, 278 K
├── neat_liq_278K-prm01-analyze.sh      # analyze script for above
├── neat_liq_278K-prm01-analyze.log     # ANALYZE9 output
├── neat_liq_298K.xyz                   # symlink → neat_liq.xyz
├── neat_liq_298K.key                   # symlink → neat_liq.key
├── neat_liq_298K.sh                    # MD run script for 298 K
├── neat_liq_298K-md.log                # MD log
├── neat_liq_298K.arc                   # full trajectory
├── neat_liq_298K.dyn                   # MD checkpoint
├── neat_liq_298K-prod.arc              # production-only arc
├── neat_liq_298K-prm01.key             # analyze key: prm perturb +Δ, 298 K
├── neat_liq_298K-prm01-analyze.sh      # analyze script for above
└── neat_liq_298K-prm01-analyze.log     # ANALYZE9 output
```

> For `N` free parameters, `2N` parameter perturbations are created per optimizer call; the analyze files are labeled `prm01` … `prm{2N}` and all submitted to the cluster simultaneously.

## Recommended Minimal Settings

For hydration free energy calculations:

| Setting | Value | Notes |
|---------|-------|-------|
| `lambda_window` | `courser` | 17 windows — good accuracy with less cost |
| `liquid_md_total_time` | 1.25 ns | Last 80% (1 ns) used in BAR analysis |
| `liquid_md_time_step` | 2.0 fs | Works well with RESPA integrator |
| `gas_md_total_time` | 1.25 ns | Last 80% (1 ns) used in BAR analysis |
| `gas_md_time_step` | 0.1 fs | Required for gas-phase stochastic dynamics |

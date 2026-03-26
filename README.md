
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
│   ├── parmOPT.py          # Parameter optimization targeting experimental HFE
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

The `utils/parmOPT.py` utility optimizes force field parameters to match experimental hydration free energies:

```bash
python utils/parmOPT.py
```

It uses `scipy.optimize.least_squares` with a central finite-difference Jacobian computed in parallel. Add these fields to `settings.yaml`:

```yaml
expt_hfe: -4.75               # Experimental HFE (kcal/mol)
opt_params: "vdwpair-401-402 3.8 0.05"  # Term, initial values
params_range: "0.5 0.02"      # Search range for each parameter
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

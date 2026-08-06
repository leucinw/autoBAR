# Phenol: HFE + dimer interaction energy (parmOPT example)

A `parmOPT.py` example that optimizes force-field parameters against **two**
targets at once, with the neat-liquid density target **disabled**:

1. **Hydration free energy (HFE)** — driven through `autoBAR.py`.
2. **Dimer interaction energy** — gas-phase `E_int = E_dimer − E_mon1 − E_mon2`
   from Tinker `analyze E`, compared against a QM reference.

Each target is optional in `parmOPT.py`; a target turns on only when its
defining key is present in `settings.yaml`. Here `expt_hfe` and `dimer_data`
are set, while `expt_density`/`expt_densities` are absent — so density is off.

## Files

| File | Purpose |
|------|---------|
| `settings.yaml` | Combined autoBAR-HFE + parmOPT settings (density omitted) |
| `amoeba09.prm` | AMOEBA parameter file |
| `phenol.xyz` | Gas-phase phenol (autoBAR `gas_xyz`) |
| `phenol_solv.xyz` | Solvated phenol box (autoBAR `box_xyz`) |
| `dimers/phenol_water.xyz` | Phenol–water H-bonded dimer for the dimer target |

## What is optimized

The van der Waals `R` and `epsilon` of the phenol hydroxyl oxygen (class 36):

```
opt_params:   "vdw-36 3.4050 0.1100"
params_range: "0.20   0.05"          # ±0.20 Å on R, ±0.05 on epsilon
```

## Run it

```bash
cd examples/Phenol-HFE-Dimer
python /path/to/autoBAR/utils/parmOPT.py
```

Progress and a per-step target-vs-current table are written to `parmOPT.log`;
the optimized parameter file is written to `amoeba09.prm.final`.

## ⚠️ Placeholder dimer data — replace before real use

`dimers/phenol_water.xyz` is a **constructed** starting geometry (O–H···O
hydrogen bond at ~1.9 Å), **not** QM-optimized, and the reference energy in
`settings.yaml`:

```yaml
dimer_data:
  - "dimers/phenol_water.xyz  -6.9  1.0"
```

uses an **illustrative** `-6.9 kcal/mol`. Both the geometry and the QM value are
placeholders to demonstrate the file format — substitute your own QM-optimized
dimer geometry and interaction energy before fitting for production.

`dimer_frag1_natoms: 13` tells parmOPT the first 13 atoms are fragment 1
(phenol) and the remaining atoms are fragment 2 (water); parmOPT splits the two
monomers automatically.

## Adding the density target

To fit neat-liquid density as a third target, add `expt_density` (with
`temperature`), a `liquid_dir` containing the neat-liquid `.xyz`, and
`production_time`. See `utils/settings.yaml` for the complete list of
density/liquid-MD keys.

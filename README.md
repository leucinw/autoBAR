
# autoBAR

An automated tool for alchemical free energy simulation using polarizable AMOEBA and AMOEBA+ force fields with Tinker software packages.

# Prerequisite
- python modules: `yaml` and `numpy`
- compiled software: `Tinker` (for gas phase simulations) and `Tinker9` (for liquid/solution phase NPT/NVT simulations)

# How to setup 

Prepare the following 4 files:
* `gas_xyz`: ligand tinker xyz file
* `box_xyz`: ligand or ligand-protein in water
* `parameters`: Tinker parameter file (see [Notes](#notes) for one-step perturbation function)
* `settings.yaml`: settings read by autoBAR.py program. Please refer to the [example file here](https://github.com/leucinw/autoBAR/blob/main/dat/settings.yaml).

# How to run 

To make it flexible to use, this program was designed to be run in either interactive or automated mode. 
In the interactive mode, an individual step can be run depending on the requirement. 
In the automated mode, this program will automatically go through all steps until it exits.

* Interactive mode
  ```shell
  # Run `setup`: generate the necessary input files for Tinker
  python autoBAR.py setup
  # Run `dynamic`: do molecular dynamics simulations at a series of lambda using Tinker/Tinker9
  python autoBAR.py dynamic
  # Run `bar`: do bar analysis using Tinker/Tinker9 after the above MD jobs finish
  python autoBAR.py bar
  # Run `result`: summarize and printout the bar analysis result
  python autoBAR.py result
  ```
* Automated mode
  ```shell
  # Run `auto`: automatically run all the above commands
  python autoBAR.py auto
  ```

# Examples

Two example systems are located in `examples` folder. They should be easy to read and understand. 

# Notes

* Suggested settings for HFE simulations
  * `lambda_window`: courser
  * `liquid_md_total_time`: 1.25 ns
  * `liquid_md_time_step`: 2.0 fs
  * `gas_md_total_time`: 1.25 ns
  * `gas_md_time_step`: 0.1 fs
  * `polar_eps`: 10e-5

* One-step perturbation is supported 
  * A `{fname}.prm_` file with small parameter perturbation need to be in the working directory
  * No need modify the `settings.yaml` file
  * This will be treated as the end state (two end states in HFE)

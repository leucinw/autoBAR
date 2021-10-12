# autoBAR
An automated tool for alchemical free energy simulation using polarizable AMOEBA and AMOEBA+ force fields with Tinker packages.

# Prerequisite
- python modules: `yaml` and `numpy`
- compiled software: `Tinker` and `Tinker9` (if GPU will be used)

# How to setup 
Prepare the following 4 files:
* `gas_xyz`: ligand tinker xyz file
* `box_xyz`: ligand or ligand-protein in water
* `tinker_prm`: Tinker parameter file
* `settings.yaml`: settings read by autoBAR.py program. [see an example file here](https://github.com/leucinw/autoBAR/blob/main/dat/settings.yaml).

# How to run 
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
# Notes
* Suggested settings for Hydration FE
  * `lambda_window`: courser
  * `liquid_md_total_time`: 1.25 ns
  * `liquid_md_time_step`: 2.0 fs
  * `gas_md_total_time`: 1.25 ns
  * `gas_md_time_step`: 0.1 fs
  * `polar_eps`: 10e-5

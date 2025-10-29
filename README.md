
# autoBAR: Free Energy Simulator with Tinker Software Packages 

A light-weighted automation tool for alchemical free energy simulations with the [Bennett Acceptance Ratio (BAR)](https://en.wikipedia.org/wiki/Bennett_acceptance_ratio) method using the polarizable AMOEBA and AMOEBA+ force fields.

# Prerequisite
- python modules: `ruamel.yaml` and `numpy`. (Hint: one can use `pip` or `conda` to install easily)
- compiled software: `Tinker` and `Tinker9`. (Note: Nothing need to do for Ren lab members).

# How to setup 

Prepare the following 4 files:
* `gas_xyz`: ligand tinker xyz file
* `box_xyz`: ligand or ligand-protein in water, with box info. in the second line
* `parameters`: Tinker parameter file (see [Notes](#notes) for one-step perturbation function)
* `settings.yaml`: settings read by autoBAR.py program. Please refer to the [example file here](https://github.com/leucinw/autoBAR/blob/main/dat/settings.yaml).

__Note: User customized tinker key files are supported (see [here](https://github.com/leucinw/autoBAR/blob/main/dat/settings.yaml))__

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
* Advanced feature
	```shell
	# Run `opt`: use one-step BAR to get numerical gradient of the FE w.r.t. tuning parameters, and optimize the FF parameters. The settings.yaml needs to be setup carefully!
	python autoBAR.py opt
	```

# Examples

Two example systems are located in `examples` folder. They should be easy to read and understand.

* For the __Ion-HFE__ system, only the `liquid` phase is necessary
	```
	├── liquid
	├── Na-water.xyz
	├── Na.xyz
	├── result.txt
	├── settings.yaml
	└── water03.prm
	```

* For the __Phenol-HFE__ system, both the `liquid` and `gas` phases are needed
	```
	├── amoeba09.prm
	├── gas
	├── liquid
	├── phenol_solv.xyz
	├── phenol.xyz
	├── result.txt
	└── settings.yaml
	```

# Notes

* Minimal settings for HFE simulations
  * `lambda_window`: courser # This reduces the number of windows without losing much accuracy
  * `liquid_md_total_time`: 1.25 ns # The last 4/5 of trajectories (1 ns) is used in BAR
  * `liquid_md_time_step`: 2.0 fs # Good with RESPA integrator
  * `gas_md_total_time`: 1.25 ns # The last 4/5 of trajectories (1 ns) is used in BAR
  * `gas_md_time_step`: 0.1 fs # Gas phase stochastic dynamics

* One-step perturbation is supported 
  * An `{fname}.prm_XX` file with small parameter perturbation need to be in the working directory
  * When there are multiple `{fname}.prm_XX` files, program will do multiple perturbations at the same time. (support upto 99 perturbations!)
  * No need modify the `settings.yaml` file
  * This will be treated as the end state (two end states if involving two legs)

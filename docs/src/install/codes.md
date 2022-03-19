# DFT Package Setup

## VASP

As described in the [ASE VASP documentation](https://wiki.fysik.dtu.dk/ase/ase/calculators/vasp.html#environment-variables), you will need to define a few environment variables. By default, you do not need to define the `ASE_VASP_COMMAND` (or `VASP_SCRIPT`) environment variable unless you plan to run QuAcc without on-the-fly error-handling using Custodian. However, you do need to define the `VASP_PP_PATH` environment variable that points to your pseudopotential library. Additionally, for the use of van der Waals functionals, you need to define the `ASE_VASP_VDW` environment variable to point to the `vdw_kernel.bindat` file distributed with VASP.

To run VASP with Custodian (the default behavior in QuAcc), you will need to define a `VASP_PARALLEL_CMD` environment variable that tells Custodian how to parallelize VASP. For instance, this might look something like `export VASP_PARALLEL_CMD="srun -N 2 --ntasks-per-node 64 -c 4 --cpu_bind=cores"`. Note, the VASP executables are not included in this environment variable.

By default, Custodian will assume that the VASP executables can be run with `vasp_std` or `vasp_gam` for standard or gamma-point calculations. If you need to use different executable names, create and modify the `vasp_custodian_settings.yaml` file found [here](https://github.com/arosen93/quacc/blob/main/quacc/defaults/custodian_settings/vasp_custodian_settings.yaml) and define the `VASP_CUSTODIAN_SETTINGS` environment variable, which should point to the newly created `vasp_custodian_settings.yaml` file.

## Gaussian

As noted in the [ASE documentation](https://wiki.fysik.dtu.dk/ase/ase/calculators/gaussian.html), ASE will look for executables named `g16`, `g09`, or `g03` (in that order). If you want to use a different executable name, you can define the `ASE_GAUSSIAN_COMMAND` environment variable accordingly.

## ORCA

As noted in the [ASE documentation](https://wiki.fysik.dtu.dk/ase/ase/calculators/orca.html), to use ORCA in parallel mode, you will need to define an environment variable named `ASE_ORCA_COMMAND` that is the full path to the ORCA executable.

## xTB

If you plan to use xTB with QuAcc, you will need to install the [xtb-python](https://github.com/grimme-lab/xtb-python) interface. This can be done via `conda install -c conda-force xtb-python`.

## Other Codes

For any other code, simply make sure to follow the instructions in [ASE's documentation](https://wiki.fysik.dtu.dk/ase/ase/calculators/calculators.html#supported-calculators) for setting up your desired calculator.

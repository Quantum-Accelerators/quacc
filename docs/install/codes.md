# Calculator Setup

!!! Tip "Start Simple"

    Just getting started? Try using the EMT or LJ recipes before worrying about setting up one of the calculators below.

Here, we outline how to ensure that quacc can run the quantum chemistry package of your choosing. You only need to follow the instructions for the code(s) you intend to use.

## DFTB+

If you plan to use DFTB+ with quacc, you will need to install the code as follows:

```bash
conda install -c conda-forge dftbplus
```

## EMT

No setup needed!

## ESPRESSO

To use quacc with ESPRESSO, you will need to define multiple environment variables. This can be done as described in the section on ["Modifying Quacc Settings"](../user/settings/settings.md). Because the Quantum Espresso has many binaries, they should be defined as such:

Actually how to do it though?

The pseudopotentials path can either be passed explicitly in the input_data keywords or set with the quacc setting. `ESPRESSO_PP_PATH`. The input_data definition takes precedence. Last resort, pw.x will always try to look into the env variable `ESPRESSO_PSEUDO` if nothing is defined.

Presets can be used to use predefined parameters, k-points, pseudopotentials, etc. In this case, the user has to pass a preset string to the function, the string must point to an existing preset (quacc.calculators.espresso.presets). Order of precedence for the parameters:

1. The `input_data` dictionary passed by the user to the function.
2. The preset parameters passed by the user if any.
3. Default parameters in some recipes.

This means that in any case the `input_data` keyword can be used to override default/preset parameters.

For other binaries, normal Namelist cards still must be passed as a dictionary to the `input_data` keyword. How to pass the additional cards will be explained in the documentation of the respective function.

```bash

## Gaussian

To use quacc with Gaussian, you will need to define the `GAUSSIAN_CMD` setting to be the path of the Gaussian executable (or the name of the executable if it is already in your `PATH`). This can be done as described in the section on ["Modifying Quacc Settings"](../user/settings/settings.md), such as by defining the following environment variable:

```bash
export QUACC_GAUSSIAN_CMD="/path/to/g16"
```

## GULP

To use quacc with GULP, you will need to download and compile GULP 6.1.2+ [per the official manual](https://gulp.curtin.edu.au/download.html). Then you will define the `GULP_CMD` setting to be the path of the GULP executable and the `GULP_LIB` setting to be the path to the GULP force field library. This can be done as described in the section on ["Modifying Quacc Settings"](../user/settings/settings.md), such as by defining the following environment variables in your `~/.bashrc`:

```bash
export QUACC_GULP_CMD="/path/to/gulp"
export QUACC_GULP_LIB="/path/to/gulp-#.#.#/Libraries"
```

!!! Tip "Receive a Compilation Error?"

    If you receive an error upon compilation, refer to [this forum post](https://matsci.org/t/installing-gulp/43158/18?u=arosen).

## Lennard Jones

No setup needed!

## NewtonNet

If you plan to use NewtonNet with Quacc, you will need to install it prior to use. This can be done as follows:

```bash
pip install quacc[newtonnet]
```

## ORCA

To use quacc with ORCA, you will need to define the `ORCA_CMD` setting to be the full, absolute path to your ORCA executable. This can be done as described in the section on ["Modifying Quacc Settings"](../user/settings/settings.md), such as by defining the following environment variable in your `~/.bashrc`:

```bash
export QUACC_ORCA_CMD="/path/to/orca/orca"
```

## Psi4

If you plan to use Psi4 with quacc, you will need to install it prior to use. This can be done as described in the [Psi4 installation guide](https://psicode.org/installs/latest/).

## Q-Chem

If you plan to use Q-Chem with Quacc, you will need to install `openbabel`. This can be done as follows:

```bash
conda install -c conda-forge openbabel
```

## tblite

If you plan to use tblite with quacc, you will need to install the tblite interface with ASE support.

```bash
pip install quacc[tblite] # only on Linux
```

## VASP

To use quacc with VASP, you will need to define several environment variables, as described in the section on ["Modifying Quacc Settings"](../user/settings/settings.md). The most important are listed below:

```bash
export QUACC_VASP_PARALLEL_CMD="srun -N 2 --ntasks-per-node 24"
export QUACC_VASP_PP_PATH="/path/to/POTCARs"
export QUACC_VASP_VDW="/path/to/directory/containing/kernel"
```

The `VASP_PARALLEL_CMD` setting tells Custodian and/or ASE how to parallelize VASP. Note that it does not include the executable.

The `VASP_PP_PATH` setting should point to the directory containing your VASP pseudopotentials. There should be two subdirectories name `potpaw_PBE` and `potpaw` for the PBE and LDA pseudopotentials, respectively. If your pseudopotential directories have a different name, create a symbolic link with the required naming scheme. We recommend setting `QUACC_VASP_PP_PATH` in your `~/.bashrc` file since this rarely changes.

The `VASP_VDW` environment variable is necessary if you are using a vdW functional and should point to the directory that contains the `vdw_kernel.bindat` file distributed with VASP. We also recommend including this in your `~/.bashrc` file since this rarely changes.

Additional settings can be specified as well, such as the name of the VASP executables if they differ from the default values (i.e. `vasp_std`, `vasp_gam`).

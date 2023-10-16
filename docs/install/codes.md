# Calculator Setup

!!! Tip

    Just getting started? Try using the EMT or LJ recipes before worrying about setting up one of the calculators below.

Here, we outline how to ensure that quacc can run the quantum chemistry package of your choosing. You only need to follow the instructions for the code(s) you intend to use.

## DFTB+

If you plan to use DFTB+ with quacc, you will need to install the code as follows:

```bash
conda install -c conda-forge dftbplus
```

## EMT

No setup needed!

## Gaussian

To use quacc with Gaussian, you will need to define the `GAUSSIAN_CMD` setting to be the path of the Gaussian executable (or the name of the executable if it is already in your `PATH`). This can be done as described in the section on ["Modifying Quacc Settings"](../user/settings/settings.md), such as by defining the following environment variable:

```bash
export QUACC_GAUSSIAN_CMD="/path/to/g16"
```

## GULP

To use quacc with GULP, you will need to define the `GULP_CMD` setting to be the path of the GULP executable (or the name of the executable if it is already in your `PATH`) and the `GULP_LIB` setting to be the path to the GULP force field library. This can be done as described in the section on ["Modifying Quacc Settings"](../user/settings/settings.md), such as by defining the following environment variables:

```bash
export QUACC_GULP_CMD="/path/to/gulp"
export QUACC_GULP_LIB="/path/to/gulp-#.#.#/Libraries"
```

## Lennard Jones

No setup needed!

## NewtonNet

If you plan to use NewtonNet with Quacc, you will need to install it prior to use. This can be done as follows:

```bash
pip install git+https://github.com/ericyuan00000/NewtonNet.git
pip install quacc[newtonnet,sella]
```

## ORCA

To use quacc with ORCA, you will need to define the `ORCA_CMD` setting to be the full, absolute path to your ORCA executable. This can be done as described in the section on ["Modifying Quacc Settings"](../user/settings/settings.md), such as by defining the following environment variable:

```bash
export QUACC_ORCA_CMD="/path/to/orca/orca"
```

## Psi4

If you plan to use Psi4 with quacc, you will need to install it prior to use. This can be done as described in the [Psi4 installation guide](https://psicode.org/installs/latest/).

## Q-Chem

If you plan to use Q-Chem with Quacc, you will need to install `openbabel` and `sella` (recommended) prior to use. This can be done as follows:

```bash
conda install -c conda-forge openbabel
pip install quacc[sella]
```

## tblite

If you plan to use tblite with quacc, you will need to install the tblite interface with ASE support.

```bash
pip install quacc[tblite] # only on Linux
```

## VASP

To use quacc with VASP, you will need to define several environment variables, as described in the section on ["Modifying Quacc Settings"](../user/settings/settings.md). The most important are listed below:

```bash
QUACC_VASP_PARALLEL_CMD="srun -N 2 --ntasks-per-node 24"
QUACC_VASP_PP_PATH="/path/to/POTCARs"
QUACC_VASP_VDW="/path/to/directory/containing/kernel"
```

The `VASP_PARALLEL_CMD` setting tells Custodian and/or ASE how to parallelize VASP. Note that it does not include the executable.

The `VASP_PP_PATH` setting should point to the directory containing your VASP pseudopotentials. There should be two subdirectories name `potpaw_PBE` and `potpaw` for the PBE and LDA pseudopotentials, respectively. If your pseudopotential directories have a different name, create a symbolic link with the required naming scheme. We recommend setting `QUACC_VASP_PP_PATH` in your `~/.bashrc` file since this rarely changes.

The `VASP_VDW` environment variable is necessary if you are using a vdW functional and should point to the directory that contains the `vdw_kernel.bindat` file distributed with VASP. We also recommend including this in your `~/.bashrc` file since this rarely changes.

Additional settings can be specified as well, such as the name of the VASP executables if they differ from the default values (i.e. `vasp_std`, `vasp_gam`).

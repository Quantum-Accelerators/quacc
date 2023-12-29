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

## MLPs

Several pre-trained "universal" machine-learned potentials (MLPs) are natively supported in quacc, including those based on [MACE](https://github.com/ACEsuit/mace), [CHGNet](https://github.com/CederGroupHub/chgnet), and [M3GNet](https://github.com/materialsvirtuallab/matgl).

To use these potentials, you will need to install the corresponding packages. This can be done as follows:

```bash
pip install quacc[ml]
pip install git+https://github.com/ACEsuit/mace.git
```

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

## Quantum ESPRESSO

To use quacc with Quantum ESPRESSO, you will first need to download and compile [Quantum ESPRESSO](https://www.quantum-espresso.org/) or install it via conda as `conda install -c conda-forge qe`. You will also need to download the relevant [pseudopotentials](https://www.quantum-espresso.org/pseudopotentials/).

Additionally, you will need to define multiple environment variables. This can be done as described in the section on ["Modifying Quacc Settings"](../user/settings/settings.md).

At minimum, you should define the `ESPRESSO_PSEUDO` setting:

```bash
QUACC_ESPRESSO_PSEUDO="/path/to/pseudopotentials"
```

Additionally, the various ESPRESSO binaries should be present in your `PATH`, or you should modify the `ESPRESSO_BIN_PATHS` quacc setting accordingly.

## TBLite

If you plan to use TBLite with quacc, you will need to install the tblite interface with ASE support.

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

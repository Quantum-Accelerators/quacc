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

As noted in the [ASE documentation](https://wiki.fysik.dtu.dk/ase/ase/calculators/gaussian.html), you will need to define an environment variable named `ASE_GAUSSIAN_COMMAND`. It should be formatted as follows:

```bash
ASE_GAUSSIAN_COMMAND="/path/to/my/gaussian_executable Gaussian.com > Gaussian.log"
```

## GULP

As noted in the [ASE documentation](https://wiki.fysik.dtu.dk/ase/ase/calculators/gulp.html), you must set the environment variables `GULP_LIB` and `ASE_GULP_COMMAND` as follows:

```bash
GULP_LIB="/path/to/my/gulp-#.#.#/Libraries"
ASE_GULP_COMMAND="/path/to/my/gulp-#.#.#/Src/gulp < gulp.gin > gulp.got"
```

## Lennard Jones

No setup needed!

## NewtonNet

If you plan to use NewtonNet with Quacc, you will need to install it prior to use. This can be done as follows:

```bash
pip install git+https://github.com/ericyuan00000/NewtonNet.git
pip install quacc[newtonnet, optimizers]
```

## ORCA

To use quacc with ORCA, you will need to define the `ORCA_CMD` setting to be the full, absolute path to your ORCA executable. This can be done as described in the section on ["Modifying Quacc Settings"](../user/settings/settings.md), such as by defining the following environment variable:

```bash
QUACC_ORCA_CMD=/path/to/my/orca/orca
```

## Psi4

If you plan to use Psi4 with quacc, you will need to install it prior to use. This can be done as described in the [Psi4 installation guide](https://psicode.org/installs/latest/).

## Q-Chem

If you plan to use Q-Chem with Quacc, you will need to install `openbabel` and `sella` (recommended) prior to use. This can be done as follows:

```bash
conda install -c conda-forge openbabel
pip install quacc[optimizers]
```

## tblite

If you plan to use tblite with quacc, you will need to install the tblite interface with ASE support.

If you are using a Linux mmachine, run the following:

```bash
pip install quacc[tblite]
```

If you are using a Mac, run the following:

```bash
conda install -c conda-forge tblite-python
```

## VASP

To use VASP with quacc, you will need to do the following, as described in greater detail in the [ASE documentation](https://wiki.fysik.dtu.dk/ase/ase/calculators/vasp.html#pseudopotentials):

- Define the `VASP_PP_PATH` environment variable that points to your pseudopotential library. We recommend including this in your `~/.bashrc` file since this rarely changes.
- If you wish to use vdW functionals, define the `ASE_VASP_VDW` environment variable to point to the `vdw_kernel.bindat` file distributed with VASP. We recommend including this in your `~/.bashrc` file since this rarely changes.

To run VASP with Custodian, you will also need to modify the relevant [quacc configuration settings](../user/settings/settings.md), most notably:

- Define the `VASP_PARALLEL_CMD` quacc setting that tells Custodian how to parallelize VASP, such as by defining an environment variable `QUACC_VASP_PARALLEL_CMD="srun -N 2 --ntasks-per-node 24"`. Note, the VASP executables are not included in this environment variable.

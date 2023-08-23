# Calculator Setup

!!! Tip

    Just getting started? Try using the EMT or LJ recipes before worrying about setting up one of the calculators below.

Here, we outline how to ensure that quacc can run the quantum chemistry package of your choosing. You only need to follow the instructions for the code(s) you intend to use.

## DFTB+

!!! Note

    [DFTB+](https://dftbplus.org/) is especially useful for periodic GFN-xTB calculations and the DFTB+ method based on Slater-Koster parameters.

If you plan to use DFTB+ with quacc, you will need to install the code as follows:

```bash
conda install -c conda-forge dftbplus
```

## EMT

!!! Note

    [Effective medium theory (EMT)](https://doi.org/10.1016/0039-6028(96)00816-3) is a semi-empirical method for modeling solids that is predominantly used for prototyping workflows. Because it is solely for demonstration purposes, it only supports the following metals: Al, Ni, Cu, Pd, Ag, Pt, and Au.

No setup needed!

## Gaussian

!!! Note

    [Gaussian](https://gaussian.com/) is an extremely popular molecular DFT code that is quite robust and easy to use.

As noted in the [ASE documentation](https://wiki.fysik.dtu.dk/ase/ase/calculators/gaussian.html), you will need to define an environment variable named `ASE_GAUSSIAN_COMMAND`. It should be formatted as follows:

```bash
ASE_GAUSSIAN_COMMAND="/path/to/my/gaussian_executable Gaussian.com > Gaussian.log"
```

## GULP

!!! Note

    [GULP](https://gulp.curtin.edu.au/) is especially useful for periodic GFN-FF calculations and force field methods. GULP can be downloaded and installed [here](https://gulp.curtin.edu.au/download.html).

As noted in the [ASE documentation](https://wiki.fysik.dtu.dk/ase/ase/calculators/gulp.html), you must set the environment variables `GULP_LIB` and `ASE_GULP_COMMAND` as follows:

```bash
GULP_LIB="/path/to/my/gulp-#.#.#/Libraries" ASE_GULP_COMMAND="/path/to/my/gulp-#.#.#/Src/gulp < gulp.gin > gulp.got"
```

## Lennard Jones

!!! Note

    [Lennard Jones (LJ)](https://en.wikipedia.org/wiki/Lennard-Jones_potential) is an empirical potential that is predominantly used for prototyping workflows for molecules.

No setup needed!

## NewtonNet

!!! Note

    [NewtonNet](https://github.com/ericyuan00000/NewtonNet) is a message passing networking for deep learning of interatomic potentials and forces, as described [here](https://pubs.rsc.org/en/content/articlehtml/2022/dd/d2dd00008c).

If you plan to use NewtonNet with Quacc, you will need to install it prior to use. This can be done via `pip install "quacc[newtonnet] @ git+https://github.com/Quantum-Accelerators/quacc.git"`.

## ORCA

!!! Note

    [ORCA](https://orcaforum.kofo.mpg.de/app.php/portal) is a free code that is especially useful for molecular DFT calculations with recently developed methods. ORCA can be downloaded and installed [here](https://orcaforum.kofo.mpg.de/app.php/dlext/).

To use quacc with ORCA, you will need to define the `ORCA_CMD` setting to be the full, absolute path to your ORCA executable. This can be done as described in the section on ["Modifying Quacc Settings"](../user/settings.md), such as by defining the following environment variable:

```bash
QUACC_ORCA_CMD=/path/to/my/orca/orca
```

## Psi4

!!! Note

    [Psi4](https://github.com/psi4/psi4) is an open-source quantum chemistry electronic structure package.

If you plan to use Psi4 with quacc, you will need to install it prior to use. This can be done as described in the [Psi4 installation guide](https://psicode.org/installs/latest/).

## tblite

!!! Note

    [tblite](https://github.com/tblite/tblite) is a code that interfaces with the xtb package for running GFN-xTB calculations.

If you plan to use tblite with quacc, you will need to install the tblite interface with ASE support.

If you are using a Linux mmachine, run the following:

```bash
pip install "quacc[tblite] @ git+https://github.com/Quantum-Accelerators/quacc.git"
```

If you are using a Mac, run the following:

```bash
conda install -c conda-forge tblite-python
```

## VASP

!!! Note

    [VASP](https://www.vasp.at/) is a very widely used code for plane-wave, periodic DFT calculations. Quacc has built-in support for automatically fixing failed VASP jobs via [Custodian](https://github.com/materialsproject/custodian).

To use VASP with quacc, you will need to do the following, as described in greater detail in the [ASE documentation](https://wiki.fysik.dtu.dk/ase/ase/calculators/vasp.html#pseudopotentials):

- Define the `VASP_PP_PATH` environment variable that points to your pseudopotential library. We recommend including this in your `~/.bashrc` file since this rarely changes.
- If you wish to use vdW functionals, define the `ASE_VASP_VDW` environment variable to point to the `vdw_kernel.bindat` file distributed with VASP. We recommend including this in your `~/.bashrc` file since this rarely changes.

To run VASP with Custodian, you will also need to modify the following [quacc configuration setting](../user/settings.md):

- Define the `VASP_PARALLEL_CMD` quacc setting that tells Custodian how to parallelize VASP, such as by defining an environment variable `QUACC_VASP_PARALLEL_CMD="srun -N 2 --ntasks-per-node 24"`. Note, the VASP executables are not included in this environment variable.

## Other Codes

For any other code, make sure to follow the instructions in [ASE's documentation](https://wiki.fysik.dtu.dk/ase/ase/calculators/calculators.html#supported-calculators) for setting up your desired calculator.

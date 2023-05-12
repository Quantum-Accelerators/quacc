# Quantum Chemistry Package Setup

Here, we outline how to ensure that QuAcc can run the quantum chemistry package of your choosing. You only need to follow the instructions for the code(s) you intend to use.

Note that some of these codes do not natively run on Windows. If you are on a Windows machine, we recommend installing and using the [Windows Subsystem for Linux (WSL)](https://ubuntu.com/wsl).

## DFTB+

DFTB+ is especially useful for periodic GFN-xTB calculations and the DFTB+ method based on Slater-Koster parameters.

If you plan to use DFTB+ with QuAcc, you will need to install the code from [here](https://dftbplus.org/). This can be done via `conda install -c conda-forge dftbplus`.

## Gaussian

Gaussian is an extremely popular molecular DFT code that is quite robust and easy to use.

As noted in the [ASE documentation](https://wiki.fysik.dtu.dk/ase/ase/calculators/gaussian.html), you will need to define an environment variable named `ASE_GAUSSIAN_COMMAND`. It should be of the form `ASE_GAUSSIAN_COMMAND="/path/to/my/gaussian_executable Gaussian.com > Gaussian.log"`.

## GULP

GULP is especially useful for periodic GFN-FF calculations and force field methods. GULP can be downloaded and installed [here](https://gulp.curtin.edu.au/).

As noted in the [ASE documentation](https://wiki.fysik.dtu.dk/ase/ase/calculators/gulp.html), you must set the environment variables `GULP_LIB="/path/to/my/gulp-#.#.#/Libraries"` and `ASE_GULP_COMMAND="/path/to/my/gulp-#.#.#/Src/gulp < gulp.gin > gulp.got"` based on where you installed GULP.

## ORCA

ORCA is a free code that is especially useful for molecular DFT calculations with recently developed methods. ORCA can be downloaded and installed [here](https://orcaforum.kofo.mpg.de/app.php/portal).

As noted in the [ASE documentation](https://wiki.fysik.dtu.dk/ase/ase/calculators/orca.html), to use ORCA in parallel mode, you will need to define an environment variable named `ASE_ORCA_COMMAND`. It should be of the form `ASE_ORCA_COMMAND="/path/to/my/orca orca.inp > orca.out"`.

## Psi4

Psi4 is especially useful for constructing and testing out new functionals, like DeepMind's DM21.

If you plan to use [Psi4](https://github.com/psi4/psi4) with QuAcc, you will need to install it prior to use. This can be done via `conda install -c psi4 psi4`.

## tblite

tblite is a code that interfaces with the xtb package for running GFN-xTB calculations. If you want to run GFN-xTB calculations on molecular systems, this is the easiest route to do so.

If you plan to use tblite with QuAcc, you will need to install the [tblite](https://github.com/tblite/tblite) interface with ASE support. This can be done via `pip install tblite[ase]`.

## VASP

VASP is a very widely used code for plane-wave, periodic DFT calculations. To use VASP with QuAcc, you will need to do the following:

1. Define the `VASP_PP_PATH` environment variable that points to your pseudopotential library. See the [ASE documentation](https://wiki.fysik.dtu.dk/ase/ase/calculators/vasp.html#pseudopotentials) for full details. We recommend including this in your `~/.bashrc` file.
2. To use vdW functionals, define the `ASE_VASP_VDW` environment variable to point to the `vdw_kernel.bindat` file distributed with VASP. Again, refer to the [ASE documentation](https://wiki.fysik.dtu.dk/ase/ase/calculators/vasp.html#pseudopotentials) for additional details. We recommend including this in your `~/.bashrc` file.
3. To run VASP with Custodian (the default behavior in QuAcc), you will need to define a `VASP_PARALLEL_CMD` environment variable that tells Custodian how to parallelize VASP. For instance, this might look something like `export VASP_PARALLEL_CMD="srun -N 2 --ntasks-per-node 24"`. Note, the VASP executables are not included in this environment variable. For convenience, we recommend specifying this environment variable at runtime so you can easily modify it.
4. By default, Custodian will assume that the VASP executables can be run with `vasp_std` or `vasp_gam` for standard or gamma-point calculations. If you need to use different executable names or wish to change any Custodian settings from the selected defaults, refer to the section on modifying QuAcc settings.

## Other Codes

For any other code, simply make sure to follow the instructions in [ASE's documentation](https://wiki.fysik.dtu.dk/ase/ase/calculators/calculators.html#supported-calculators) for setting up your desired calculator.

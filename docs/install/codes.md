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

To use quacc with GULP, you will need to download and compile GULP 6.2+ [per the official manual](https://gulp.curtin.edu.au/download.html). Then you will define the `GULP_CMD` setting to be the path of the GULP executable and the `GULP_LIB` setting to be the path to the GULP force field library. This can be done as described in the section on ["Modifying Quacc Settings"](../user/settings/settings.md), such as by defining the following environment variables in your `~/.bashrc`:

```bash
export QUACC_GULP_CMD="/path/to/gulp"
export QUACC_GULP_LIB="/path/to/gulp-#.#.#/Libraries"
```

??? Tip "Receive a Compilation Error?"

    If you receive an error upon compilation, you likely are using an old version of gfortran. Try `./mkgulp_old_gfortran` instead of `./mkgulp` in the `Src` directory.

## Lennard Jones

No setup needed!

## Foundation MLIPs

Several pre-trained "foundation" machine-learned interatomic potentials (MLIPs) are natively supported in quacc. Quacc supports all MLIPs provided by [rootstock](https://github.com/Garden-AI/rootstock), [matcalc](https://github.com/materialyzeai/matcalc/tree/main#supported-foundation-potentials), and [FAIR-Chem](https://github.com/FAIR-Chem/fairchem).

To use these potentials, you will need to install the corresponding packages:

```bash
pip install quacc[mlip]
```

## ONETEP

If you plan on using ONETEP with quacc, you will need to obtain and install ONETEP. This can be done as described in the [How to Get Onetep](https://onetep.org/code/) section of the software documentation.

At minimum, you will need to define and set `ONETEP_CMD` to be the full path to your ONETEP binary. You can also specify the `ONETEP_PP_PATH` to be the path to the pseudopotentials. This can be done as described in the section on ["Modifying Quacc Settings"](../user/settings/settings.md).

An example is provided below on how to define the commands in your `~/.bashrc`:

```bash
export ONETEP_CMD="/path/to/onetep/binary"
export ONETEP_PP_PATH="/path/to/my/pseudos"
```

## ORCA

To use quacc with ORCA, you will need to define the `ORCA_CMD` setting to be the full, absolute path to your ORCA executable. This can be done as described in the section on ["Modifying Quacc Settings"](../user/settings/settings.md), such as by defining the following environment variable in your `~/.bashrc`:

```bash
export QUACC_ORCA_CMD="/path/to/orca/orca"
```

## Psi4

If you plan to use Psi4 with quacc, you will need to install it prior to use. This can be done as follows:

```bash
conda install -n base conda-libmamba-solver
conda install psi4 -c conda-forge --solver libmamba
```

## Q-Chem

If you plan to use Q-Chem with Quacc, you will need to install `openbabel`. This can be done as follows:

```bash
conda install -c conda-forge openbabel
```

## Quantum ESPRESSO

To use quacc with Quantum ESPRESSO, you will first need to download and install [Quantum ESPRESSO](https://www.quantum-espresso.org/). This can be most easily down as follows:

```bash
conda install -c conda-forge qe
```

You will also need to download the relevant [pseudopotentials](https://www.materialscloud.org/discover/sssp/table/efficiency).

Finally, you will need to define multiple environment variables. This can be done as described in the section on ["Modifying Quacc Settings"](../user/settings/settings.md).

At minimum, you should define the `ESPRESSO_PSEUDO` setting:

```bash
export QUACC_ESPRESSO_PSEUDO="/path/to/pseudopotentials"
```

The various ESPRESSO binaries should be present in your `PATH`, or you should modify the `ESPRESSO_BIN_PATHS` quacc setting accordingly.

## TBLite

If you plan to use TBLite with quacc, you will need to install the tblite interface with ASE support.

=== "Conda"

    ```bash
    conda install -c conda-forge tblite-python
    ```

=== "Pip"

    ```bash
    pip install quacc[tblite] # only on Linux
    ```

Refer to the [TBLite documentation](https://tblite.readthedocs.io/en/latest/tutorial/parallel.html) for details on how to parallelize calculations and adjust the memory limits.

## VASP

To use quacc with VASP, you will need to define several environment variables, as described in the section on ["Modifying Quacc Settings"](../user/settings/settings.md). The most important are listed below:

```bash
export QUACC_VASP_PARALLEL_CMD="srun -N 2 --ntasks-per-node 24"
export QUACC_VASP_PP_PATH="/path/to/POTCARs"
export QUACC_VASP_PP_VERSION="64"
```

The `VASP_PARALLEL_CMD` setting tells Custodian and/or ASE how to parallelize VASP. Note that it does not include the executable.

The `VASP_PP_PATH` setting should be defined as described in the [ASE documentation](https://wiki.fysik.dtu.dk/ase/ase/calculators/vasp.html#pseudopotentials) such that it is a directory containing each VASP pseudopotential folder. Similarly, the `VASP_PP_VERSION` setting defines the default version of the pseudopotentials to use if none is specified in a chosen recipe.

Additional settings can be specified as well, such as the name of the VASP executables if they differ from the default values (i.e. `vasp_std`, `vasp_gam`).

### Combined SCF and NSCF calculations

With ASE's [KPOINTS_OPT support](https://gitlab.com/ase/ase/-/merge_requests/4298)
and VASP 6.3 or newer, `non_scf_job` can perform both steps in one invocation:

```python
from ase.build import bulk
from quacc.recipes.vasp.core import non_scf_job

result = non_scf_job(
    bulk("Si", "diamond", a=5.43),
    use_kpoints_opt=True,
    kpts_mode="line",
    kpts=(6, 6, 6),
)
```

The preset and `kpts`, `kspacing`, or `pmg_kpts` configure the uniform SCF grid.
The NSCF grid is generated from `kpts_mode`, `uniform_kppvol`, or
`line_kpt_density`. Override it with `kpts_opt` and, for explicit reciprocal
coordinates, `reciprocal_opt=True`. With VASP 6.6 or newer, pass
`kspacing_opt=0.1` instead to generate the optional mesh in VASP.
The standard VASP executable is used for optional grids, including when the
SCF grid contains only Gamma.

No previous calculation is required in this mode. If `prev_dir` is supplied,
its charge density and wavefunctions are copied and `nbands_factor` scales
the prior NBANDS. Set `nbands` explicitly when starting fresh if more empty
bands are needed. For unscreened hybrid functionals, use `hfrcut=-1` as
described in the [VASP documentation](https://vasp.at/wiki/KPOINTS_OPT).
ASE calculator results describe the SCF step; optional electronic results
depend on the output parser. This feature currently requires the upstream
ASE change linked above and has not been validated with a real VASP run.

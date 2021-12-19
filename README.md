# HT-ASE
Various scripts enhancing ASE for high-throughput DFT.

## Installation
Install HT-ASE via `pip install .` in the base directory. We recommend doing so in a clean virtual (e.g. [Miniconda](https://docs.conda.io/en/latest/miniconda.html)) environment.

General:
- Add `~/htase` to your `$PYTHONPATH` environment variable and add
`~/htase/bin` to `$PATH` (assuming `~/htase` is where your HT-ASE folder is).

For VASP support:
- Set the `VASP_PP_PATH` environment variable to point to your library of VASP PAW pseudopotentials, as described in the [ASE documentation](https://wiki.fysik.dtu.dk/ase/ase/calculators/vasp.html#pseudopotentials).
- Set the `ASE_VASP_VDW` environment variable to point to your VASP vdW kernel file (typically named `vdw_kernel.bindat`), as described in the [ASE documentation](https://wiki.fysik.dtu.dk/ase/ase/calculators/vasp.html#pseudopotentials). This is needed if vdW functionals are used.
- Set the `RUN_VASP_CUSTODIAN` environment variable to 
- Set the 

For database support:
- Set the `JOBFLOW_CONFIG_FILE` environment variable to point to your `jobflow.yaml` file, which contains information about where to store calculation outputs, as described [here](https://materialsproject.github.io/atomate2/user/install.html#jobflow-yaml). If not set, outputs will be stored in memory via the `MemoryStore`.

## Requirements
Can be installed via `pip install -r requirements.txt`:
- [ASE](https://gitlab.com/ase/ase)
- [Pymatgen](https://github.com/materialsproject/pymatgen)
- [Custodian](https://github.com/materialsproject/custodian)
- [Jobflow](https://github.com/materialsproject/jobflow)
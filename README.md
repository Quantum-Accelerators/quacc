# HT-ASE
Various scripts enhancing ASE for high-throughput DFT.

## Installation
Install HT-ASE via `pip install .` in the base directory. We recommend doing so in a clean virtual (e.g. [Miniconda](https://docs.conda.io/en/latest/miniconda.html)) environment.

In addition, you will want to define several environment variables (typically in your `~/.bashrc`) to tell HT-ASE and the underlying ASE code how to run various codes.

Required for VASP:
- Set the `VASP_PP_PATH` environment variable to point to your library of VASP PAW pseudopotentials, as described in the [ASE documentation](https://wiki.fysik.dtu.dk/ase/ase/calculators/vasp.html#pseudopotentials).
- Set the `HTASE_DIR` environment variable to point to the base directory of the HT-ASE package. This does nothing on its own. It is simply for convenience since we will reference it a lot.
- `export ASE_VASP_COMMAND="python ${HTASE_DIR}/htase/custodian/run_vasp_custodian.py"`. This tell ASE to run Custodian-powered VASP.
- Edit the `vasp_cmd` and `vasp_gamma_cmd` in the `${HTASE_DIR}/htase/custodian/vasp_custodian_settings.yaml` [file](https://github.com/arosen93/HT-ASE/blob/main/htase/custodian/vasp_custodian_settings.yaml) to tell Custodian how to run VASP on your supercomputer. The file also contains some defualt settings for Custodian. If you want different settings for various projects (e.g. different numbers of nodes, different Custodian handlers), you can make your own and define the path in an `VASP_CUSTODIAN_SETTINGS` environment variable at runtime.

Optional for VASP:
- Set the `ASE_VASP_VDW` environment variable to point to your VASP vdW kernel file (typically named `vdw_kernel.bindat`), as described in the [ASE documentation](https://wiki.fysik.dtu.dk/ase/ase/calculators/vasp.html#pseudopotentials). This is needed if vdW functionals are used.
- `export ASE_VASP_SETUPS="${HTASE_DIR}/defaults/user_setups/vasp"`. This allows you to easily access our [custom setups](https://github.com/arosen93/HT-ASE/blob/main/htase/defaults/user_setups/vasp) (e.g. `setups='$pbe54'`) when instantiating your calculator.

Required for database support:
- Make a `jobflow.yaml` as described in the [Atomate2 documentation](https://materialsproject.github.io/atomate2/user/install.html#jobflow-yaml) and then set the `JOBFLOW_CONFIG_FILE` environment variable to point to this `jobflow.yaml` file. The `jobflow.yaml` contains information about where to store calculation outputs. If the config file is not found by jobflow, serialized outputs will be stored in memory.

## Requirements
Can be installed via `pip install -r requirements.txt`:
- [ASE](https://gitlab.com/ase/ase)
- [Pymatgen](https://github.com/materialsproject/pymatgen)
- [Custodian](https://github.com/materialsproject/custodian)
- [Jobflow](https://github.com/materialsproject/jobflow)

# HT-ASE
HT-ASE enhances [ASE](https://wiki.fysik.dtu.dk/ase/index.html) for high-throughput DFT. Some features include:
- Support for running VASP in ASE via [Custodian](https://github.com/materialsproject/custodian) for on-the-fly error handling.
- A smarter ASE-based VASP calculator with a "co-pilot" mode that will automatically adjust INCAR flags that go against what is in the VASP manual. This is inspired by Pymatgen's [handling of input sets](https://github.com/materialsproject/pymatgen/blob/master/pymatgen/io/vasp/sets.py).
- Support for Pymatgen [automatic k-point generation schemes](https://pymatgen.org/pymatgen.io.vasp.inputs.html) in the ASE calculator itself.
- Easy integration with [Jobflow](https://materialsproject.github.io/jobflow/) for the simple construction of complex workflows and ability to store results in database format. By extension, this also makes it possible to easily use ASE with [Fireworks](https://github.com/materialsproject/fireworks) for job management.

In practice, the goal here is to enable the development of [Atomate2](https://github.com/materialsproject/atomate2)-like workflows centered around ASE with a focus on rapid workflow construction and prototyping. The speed of workflow development comes into play here because ASE is largely calculator-agnostic, making it possible to construct and link together workflows for dozens of simulation packages without breaking a sweat. Additionally, because ASE is mostly calculator-agnostic, rapid prototyping for new workflows can be done with semi-empirical methods (e.g. effective medium theory) before switching over to your production code of choice.

## Installation
Install HT-ASE via `pip install .` in the base directory. We recommend doing so in a clean virtual (e.g. [Miniconda](https://docs.conda.io/en/latest/miniconda.html)) environment.

In addition, you will want to define several environment variables (typically in your `~/.bashrc`), as outlined below.

**Required for VASP**:
- Set the `VASP_PP_PATH` environment variable to point to your library of VASP PAW pseudopotentials, as described in the [ASE documentation](https://wiki.fysik.dtu.dk/ase/ase/calculators/vasp.html#pseudopotentials).
- Set the `HTASE_DIR` environment variable to point to the base directory of the HT-ASE package. This does nothing on its own. It is simply for convenience since we will reference it a lot.
- `export ASE_VASP_COMMAND="python ${HTASE_DIR}/htase/custodian/run_vasp_custodian.py"`. This tell ASE to run Custodian-powered VASP.
- Edit the `vasp_cmd` and `vasp_gamma_cmd` in the `${HTASE_DIR}/htase/custodian/vasp_custodian_settings.yaml` [file](https://github.com/arosen93/HT-ASE/blob/main/htase/custodian/vasp_custodian_settings.yaml) to tell Custodian how to run VASP on your supercomputer. The file also contains some defualt settings for Custodian. If you want different settings for various projects (e.g. different numbers of nodes, different Custodian handlers), you can make your own `vasp_custodian_settings.yaml` file and define the path to it in an `VASP_CUSTODIAN_SETTINGS` environment variable at runtime.

**Optional for VASP**:
- Set the `ASE_VASP_VDW` environment variable to point to your VASP vdW kernel file (typically named `vdw_kernel.bindat`), as described in the [ASE documentation](https://wiki.fysik.dtu.dk/ase/ase/calculators/vasp.html#pseudopotentials). This is needed if vdW functionals are used.
- `export ASE_VASP_SETUPS="${HTASE_DIR}/defaults/user_setups/vasp"`. This allows you to easily access our [custom setups](https://github.com/arosen93/HT-ASE/blob/main/htase/defaults/user_setups/vasp) (e.g. `setups='$pbe54'`) when instantiating your calculator.

**Required for database support**:
- Make a `jobflow.yaml` as described in the [Atomate2 documentation](https://materialsproject.github.io/atomate2/user/install.html#jobflow-yaml) and then set the `JOBFLOW_CONFIG_FILE` environment variable to point to this `jobflow.yaml` file. The `jobflow.yaml` contains information about where to store calculation outputs. If the config file is not found by jobflow, serialized outputs will be stored in memory.

## Minimal Examples
### SmartVasp Calculator
To use HT-ASE's `SmartVasp()` calculator, simply import it from `htase.calculators.vasp` and use it with any of the [input arguments](https://wiki.fysik.dtu.dk/ase/ase/calculators/vasp.html) in a typical ASE `Vasp()` calculator. The only differences for the user are that that the first argument must be the ASE `Atoms` object, and the calculator should not be attached to the `Atoms` object since this will be done automatically.

```python
from htase.calculators.vasp import SmartVasp
from ase.build import bulk

atoms = bulk("Cu")
SmartVasp(atoms, xc="PBE", setups="$pbe54", kpts={"reciprocal_density":50})
atoms.get_potential_energy()
```
In the above example, there are already several enhancements compared to standard ASE. First, `SmartVasp()` will recognize that the requested `kpts` is an [automatic k-point generation scheme](https://pymatgen.org/pymatgen.io.vasp.inputs.html#pymatgen.io.vasp.inputs.Kpoints.automatic_density_by_vol) from Pymatgen and will use that to generate the k-points for you. `SmartVasp()` will also recognize that you have an element with d electrons and will set `LMAXMIX=4` per the VASP manual. Assuming `ASE_VASP_SETUPS` was set as outliend in the instructions, `setups="$pbe54"` will read in our [pre-defined setups library](https://github.com/arosen93/HT-ASE/tree/main/htase/defaults/user_setups/vasp) packaged with HT-ASE. Finally, when the potential energy is requested, HT-ASE will run VASP using Custodian. 

### Jobflow Integration

To convert the jobflow job/flow to a Fireworks firework/workflow, see the `jobflow.managers.fireworks.job_to_firework()` and `jobflow.managers.fireworks.flow_to_workflow()` functions in the [Jobflow documentation](https://materialsproject.github.io/jobflow/jobflow.managers.html#module-jobflow.managers.fireworks).


## Requirements
Can be installed via `pip install -r requirements.txt`:
- [ASE](https://gitlab.com/ase/ase)
- [Pymatgen](https://github.com/materialsproject/pymatgen)
- [Custodian](https://github.com/materialsproject/custodian)
- [Jobflow](https://github.com/materialsproject/jobflow)

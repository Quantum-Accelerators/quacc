# HT-ASE
HT-ASE enhances [ASE](https://wiki.fysik.dtu.dk/ase/index.html) for high-throughput DFT. Some features include:
- Support for running VASP in ASE via [Custodian](https://github.com/materialsproject/custodian) for on-the-fly error handling.
- A smarter ASE-based VASP calculator with an optional "co-pilot" mode that will automatically adjust INCAR flags if they go against what is in the [VASP manual](https://www.vasp.at/wiki/index.php/Main_page). This is inspired by Pymatgen's [handling of input sets](https://github.com/materialsproject/pymatgen/blob/master/pymatgen/io/vasp/sets.py).
- Support for Pymatgen [automatic k-point generation schemes](https://www.vasp.at/wiki/index.php/Main_page) in the ASE calculator itself.
- Easy integration with [Jobflow](https://materialsproject.github.io/jobflow/) for the simple construction of complex workflows and ability to store results in database format. By extension, this also makes it possible to easily use ASE with [Fireworks](https://github.com/materialsproject/fireworks) for job management.

In practice, the goal here is to enable the development of [Atomate2](https://github.com/materialsproject/atomate2)-like workflows centered around ASE with a focus on rapid workflow construction and prototyping. The speed of workflow development comes into play here because ASE is largely calculator-agnostic, making it possible to construct and link together workflows for dozens of simulation packages without breaking a sweat. Additionally, rapid prototyping for new workflows can be done with semi-empirical methods (e.g. effective medium theory) before switching over to your production code of choice.
<p align="center">
<img src="https://imgs.xkcd.com/comics/standards_2x.png" alt="xkcd Comic" width="528" height="300">
<p align="center">
Credit: xkcd
</p>

## Installation
0. Make sure you have Python 3.7+ installed, preferable in a clean virtual (e.g. [Miniconda](https://docs.conda.io/en/latest/miniconda.html)) environment.
1. Run the following command in a convenient place (e.g. `~/software`) to install HT-ASE:
```bash
git clone https://github.com/arosen93/htase.git && cd htase && pip install -r requirements.txt && pip install -e .
```
2. In addition, you will want to define several environment variables (e.g. in yout `~/.bashrc`), as outlined below.
```bash
export VASP_PP_PATH="/path/to/pseudopotential/library" # tells ASE where the VASP PAW pseudopotentials are
export HTASE_DIR="/path/to/htase" # path to this package (only used for convenience below)
export VASP_CUSTODIAN_SETTINGS="${HTASE_DIR}/htase/custodian/vasp_custodian_settings.yaml" # path to Custodian settings
export ASE_VASP_COMMAND="python ${HTASE_DIR}/htase/custodian/run_vasp_custodian.py" # tells ASE to run Custodian-powered VASP
export ASE_VASP_SETUPS="${HTASE_DIR}/htase/defaults/user_setups/vasp" # to access HT-ASE pseudopotential defaults (optional)
export ASE_VASP_VDW="/path/to/vdw_kernel.bindat" # for vdW functionals (optional)
export JOBFLOW_CONFIG_FILE=/path/to/jobflow.yaml # for jobflow Store support (optional). 
```

For guidance with setting up `VASP_PP_PATH` and `ASE_VASP_VDW`, see the [ASE Vasp calculator docs](https://wiki.fysik.dtu.dk/ase/ase/calculators/vasp.html#pseudopotentials). For guidance with setting up `JOBFLOW_CONFIG_FILE`, see the [Jobflow API docs](https://materialsproject.github.io/jobflow/jobflow.settings.html?highlight=jobflow_config_file#jobflow.settings.JobflowSettings).
3. Edit the `vasp_cmd` and `vasp_gamma_cmd` in the `vasp_custodian_settings.yaml` [file](https://github.com/arosen93/HT-ASE/blob/main/htase/custodian/vasp_custodian_settings.yaml) to tell Custodian how to run VASP on your supercomputer. The file also contains some defualt settings for Custodian. If you want different settings for various projects (e.g. different numbers of nodes, different Custodian handlers), you can make a new `vasp_custodian_settings.yaml` file and define the path to it in a new `VASP_CUSTODIAN_SETTINGS` environment variable at runtime.

**Required for database support**:
- Make a `jobflow.yaml` as described in the [Atomate2 documentation](https://materialsproject.github.io/atomate2/user/install.html#jobflow-yaml) and then add the following to yout `~/.bashrc`:
```bash
```
The `jobflow.yaml` contains information about where to store calculation outputs. If the config file is not found by jobflow, serialized outputs will be stored in memory.

## Minimal Examples
### SmartVasp Calculator
To use HT-ASE's `SmartVasp()` calculator, simply import it from `htase.calculators.vasp` and use it with any of the [input arguments](https://wiki.fysik.dtu.dk/ase/ase/calculators/vasp.html) in a typical ASE `Vasp()` calculator. The only differences for the user are that the first argument must be the ASE `Atoms` object, and it returns an `Atoms` object with an enhanced `Vasp()` calculator already attached.

```python
from htase.calculators.vasp import SmartVasp
from ase.build import bulk

atoms = bulk("Cu")
atoms = SmartVasp(atoms, xc="PBE", setups="$pbe54", kpts={"reciprocal_density":100})
atoms.get_potential_energy()
```
In the above example, there are already several enhancements compared to standard ASE. First, `SmartVasp()` will recognize that the requested `kpts` is an [automatic k-point generation scheme](https://pymatgen.org/pymatgen.io.vasp.inputs.html#pymatgen.io.vasp.inputs.Kpoints.automatic_density_by_vol) from Pymatgen and will use that to generate the k-points for you. `SmartVasp()` will also recognize that you have an element with d electrons and will set `LMAXMIX=4` per the VASP manual. Assuming `ASE_VASP_SETUPS` was set as outliend in the instructions, `setups="$pbe54"` will read from our [pre-defined setups library](https://github.com/arosen93/HT-ASE/tree/main/htase/defaults/user_setups/vasp) packaged with HT-ASE. Finally, when the potential energy is requested, HT-ASE will run VASP using Custodian. 

### Jobflow Integration

### Fireworks Integration
To convert the jobflow job/flow to a Fireworks firework/workflow, refer to the [Jobflow documentation](https://materialsproject.github.io/jobflow/jobflow.managers.html#module-jobflow.managers.fireworks).

## Requirements
Python 3.7+ is required in addition to the following packages:
- [ASE](https://gitlab.com/ase/ase)
- [Pymatgen](https://github.com/materialsproject/pymatgen)
- [Custodian](https://github.com/materialsproject/custodian)
- [Jobflow](https://github.com/materialsproject/jobflow)

![example workflow](https://github.com/arosen93/htase/actions/workflows/workflow.yaml/badge.svg)
[![CodeFactor](https://www.codefactor.io/repository/github/arosen93/htase/badge)](https://www.codefactor.io/repository/github/arosen93/htase)

# HT-ASE (🚧 Under Construction 🚧)
HT-ASE enhances [ASE](https://wiki.fysik.dtu.dk/ase/index.html) for high-throughput DFT. Some features include:
- Support for running VASP in ASE via [Custodian](https://github.com/materialsproject/custodian) for on-the-fly error handling.
- A smarter ASE-based VASP calculator with an optional "co-pilot" mode that will automatically adjust INCAR flags if they go against what is in the [VASP manual](https://www.vasp.at/wiki/index.php/Main_page).
- Support for Pymatgen's [automatic k-point generation schemes](https://pymatgen.org/pymatgen.io.vasp.inputs.html?highlight=kpoints#pymatgen.io.vasp.inputs.Kpoints) in the ASE calculator itself.
- The ability to read in pre-defined ASE calculators with settings defined in YAML format.
- Easy integration with [Jobflow](https://materialsproject.github.io/jobflow/) for the simple construction of complex workflows and ability to store results in database format. By extension, this also makes it possible to easily use ASE with [Fireworks](https://github.com/materialsproject/fireworks) for job management.

In practice, the goal here is to enable the development of [Atomate2](https://github.com/materialsproject/atomate2)-like workflows centered around ASE with a focus on rapid workflow construction and prototyping. The speed of workflow development comes into play because ASE is largely calculator-agnostic, making it possible to construct and link together workflows for dozens of simulation packages without breaking a sweat. Additionally, rapid prototyping for new workflows can be done with semi-empirical methods (e.g. effective medium theory) before switching over to your production code of choice.
<p align="center">
<img src="https://imgs.xkcd.com/comics/standards_2x.png" alt="xkcd Comic" width="528" height="300">
<p align="center">
Credit: xkcd
</p>

## Minimal Examples
### SmartVasp Calculator
To use HT-ASE's `SmartVasp()` calculator, simply import it from `htase.calculators.vasp` and use it with any of the [input arguments](https://wiki.fysik.dtu.dk/ase/ase/calculators/vasp.html#ase.calculators.vasp.Vasp) in a typical ASE `Vasp()` calculator. The only differences for the user are that the first argument must be the ASE `Atoms` object, and it returns an `Atoms` object with an enhanced `Vasp()` calculator already attached. There are also some newly introduced parameters, like `auto_kpts` for Pymatgen-generated k-point grids.

```python
from htase.calculators.vasp import SmartVasp
from ase.build import bulk

atoms = bulk("Cu") # example Atoms object
atoms = SmartVasp(atoms, preset="BulkRelaxSet", encut=600, auto_kpts={"reciprocal_density": 200}) # set calculator
atoms.get_potential_energy() # run VASP
```

### Jobflow Integration
```python
from htase.calculators.vasp import SmartVasp
from htase.schemas.vasp import summarize
from ase.io.jsonio import encode, decode
from ase.build import bulk
from jobflow import job, Flow
from jobflow.managers.local import run_locally

# Constrct an Atoms object
atoms = bulk("Cu") 

@job
def run_relax_static(atoms_json, static=False):

    # Decode JSON to Atoms
    atoms = decode(atoms_json)
    
    # Set calculator
    atoms = SmartVasp(atoms, preset="BulkRelaxSet")
    if static:
        atoms.calc.set(nsw=0)
        
    # Run VASP
    atoms.get_potential_energy()
    
    # Return serialized results
    return {"atoms": encode(atoms), "results": summarize.get_results()}

# Define a relax and then static job in a flow
job1 = run_relax_static(encode(atoms))
job2 = run_relax_static(job1.output["atoms"], static=True)
flow = Flow([job1, job2], output=job2.output)

# Run locally
responses = run_locally(flow, create_folders=True)
```
### Fireworks Integration
To convert the jobflow job/flow to a Fireworks firework/workflow, refer to the [Jobflow documentation](https://materialsproject.github.io/jobflow/jobflow.managers.html#module-jobflow.managers.fireworks). It can be done via
```python
from jobflow.managers.fireworks import flow_to_workflow

# Convert a flow (like the one above) to a fireworks WorkFlow object
wf = flow_to_workflow(flow)

# Submit the workflow to the FireWorks launchpad
# This is done instead of run_locally()
lpad = LaunchPad.auto_load()
lpad.add_wf(wf)
```

## Installation
0. Make sure you have Python 3.7+ installed, preferable in a clean virtual (e.g. [Miniconda](https://docs.conda.io/en/latest/miniconda.html)) environment.
1. Run the following command in a convenient place (e.g. `~/software`) to install HT-ASE:
```bash
git clone https://github.com/arosen93/htase.git && cd htase && pip install -r requirements.txt && pip install -e .
```
2. You will want to define several environment variables (e.g. in your `~/.bashrc`), as outlined below:
```bash
export VASP_PP_PATH="/path/to/pseudopotential/library" # tells ASE where the VASP PAW pseudopotentials are
export HTASE_DIR="/path/to/htase" # path to this package (only used for convenience below)
export VASP_CUSTODIAN_SETTINGS="${HTASE_DIR}/htase/custodian/vasp/vasp_custodian_settings.yaml" # path to Custodian settings
export ASE_VASP_COMMAND="python ${HTASE_DIR}/htase/custodian/vasp/run_vasp_custodian.py" # tells ASE to run Custodian-powered VASP
export ASE_VASP_VDW="/path/to/vdw_kernel.bindat" # for vdW functionals (optional)
export JOBFLOW_CONFIG_FILE=/path/to/jobflow.yaml # for jobflow Store support (optional). 
```

For guidance with setting up `VASP_PP_PATH` and `ASE_VASP_VDW`, see the [ASE Vasp calculator docs](https://wiki.fysik.dtu.dk/ase/ase/calculators/vasp.html#pseudopotentials). For guidance with setting up `JOBFLOW_CONFIG_FILE`, see the [Jobflow API docs](https://materialsproject.github.io/jobflow/jobflow.settings.html?highlight=jobflow_config_file#jobflow.settings.JobflowSettings).

3. Edit the `vasp_cmd` and `vasp_gamma_cmd` in the `vasp_custodian_settings.yaml` [file](https://github.com/arosen93/HT-ASE/blob/main/htase/custodian/vasp/vasp_custodian_settings.yaml) to tell Custodian how to run VASP on your supercomputer. The file also contains some defualt settings for Custodian. If you want different settings for various projects (e.g. different numbers of nodes, different Custodian handlers), you can make a new `vasp_custodian_settings.yaml` file and define the path to it in the `VASP_CUSTODIAN_SETTINGS` environment variable at runtime.

## Requirements
Python 3.7+ is required in addition to the following packages:
- [ASE](https://gitlab.com/ase/ase)
- [Pymatgen](https://github.com/materialsproject/pymatgen)
- [Custodian](https://github.com/materialsproject/custodian)
- [Jobflow](https://github.com/materialsproject/jobflow)
- [Atomate2](https://github.com/materialsproject/atomate2)

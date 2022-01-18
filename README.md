![example workflow](https://github.com/arosen93/htase/actions/workflows/workflow.yaml/badge.svg)
[![CodeFactor](https://www.codefactor.io/repository/github/arosen93/htase/badge)](https://www.codefactor.io/repository/github/arosen93/htase)

# HT-ASE (🚧 Under Construction 🚧)
HT-ASE enhances [ASE](https://wiki.fysik.dtu.dk/ase/index.html) for high-throughput DFT. Some features include:
- Support for running VASP in ASE via [Custodian](https://github.com/materialsproject/custodian) for on-the-fly error handling.
- A smarter ASE-based VASP calculator with an optional "co-pilot" mode that will automatically adjust INCAR flags if they go against what is in the [VASP manual](https://www.vasp.at/wiki/index.php/Main_page).
- Support for Pymatgen's [automatic k-point generation schemes](https://pymatgen.org/pymatgen.io.vasp.inputs.html?highlight=kpoints#pymatgen.io.vasp.inputs.Kpoints) in the ASE calculator itself.
- The ability to read in pre-defined ASE calculators with settings defined in YAML format.
- Easy integration with [Jobflow](https://materialsproject.github.io/jobflow/) for the simple construction of complex workflows and ability to store results in database format. By extension, this also makes it possible to easily use ASE with [Fireworks](https://github.com/materialsproject/fireworks) for job management.

In practice, the goal here is to enable the development of [Atomate2](https://github.com/materialsproject/atomate2)-like workflows centered around ASE with a focus on rapid workflow construction and prototyping. The speed of workflow development comes into play because ASE is largely calculator-agnostic, making it possible to construct and link together workflows for many different simulation packages without breaking a sweat. Additionally, rapid prototyping for new workflows can be done with semi-empirical methods (e.g. effective medium theory) before switching over to your production code of choice.
<p align="center">
<img src="https://imgs.xkcd.com/comics/standards_2x.png" alt="xkcd Comic" width="528" height="300">
<p align="center">
Credit: xkcd
</p>

## Minimal Examples
### SmartVasp Calculator
In direct analogy to the conventional way of running ASE, HT-ASE has a calculator called `SmartVasp()` that takes any of the [input arguments](https://wiki.fysik.dtu.dk/ase/ase/calculators/vasp.html#ase.calculators.vasp.Vasp) in a typical ASE `Vasp()` calculator but supports several additional keyword arguments to supercharge your workflow. It can also adjust your settings on-the-fly if they go against the VASP manual. The main differences for the seasoned ASE user are that the first argument must be an ASE `Atoms` object, and it returns an `Atoms` object with an enhanced `Vasp()` calculator already attached.

The example below runs a relaxation of bulk Cu using the RPBE functional with the remaining settings taken from a pre-defined set ("preset") of calculator input arguments.

```python
from htase.calculators.vasp import SmartVasp
from ase.build import bulk

atoms = bulk("Cu") # example Atoms object
atoms = SmartVasp(atoms, xc='rpbe', preset="BulkRelaxSet") # set calculator
atoms.get_potential_energy() # run VASP w/ Custodian
```

### Jobflow Integration
The above example can be converted to a format suitable for constructing a Jobflow flow simply by defining it in a function with a `@job` wrapper immediately preceeding it. One nuance of Jobflow is that the inputs and outputs must be JSON serializable (so that it can be easily stored in a database), but otherwise it works the same.

```python
from htase.calculators.vasp import SmartVasp
from htase.schemas.vasp import summarize
from ase.io.jsonio import decode

#-----Jobflow Function-----
@job
def run_relax(atoms_json):

    # Decode JSON to Atoms
    atoms = decode(atoms_json)
            
    # Run VASP
    atoms = SmartVasp(atoms, xc='rpbe', preset="BulkRelaxSet")
    atoms.get_potential_energy()
    
    # Return serialized results
    results = summarize.get_results(atoms)
    return results
```
```python
from ase.build import bulk
from ase.io.jsonio import encode
from jobflow import job, Flow
from jobflow.managers.local import run_locally

#-----Make and Run a Flow-----
# Constrct an Atoms object
atoms = bulk("Cu")

# Call the relaxation function
job1 = run_relax(encode(atoms))

# Define the flow
flow = Flow([job1])

# Run locally
responses = run_locally(flow, create_folders=True)
```
### Fireworks Integration
For additional details on how to convert a Jobflow job or flow to a Fireworks firework or workflow, refer to the [Jobflow documentation](https://materialsproject.github.io/jobflow/jobflow.managers.html#module-jobflow.managers.fireworks). 

## Installation
1. Make sure you have Python 3.7+ installed, preferable in a clean virtual (e.g. [Miniconda](https://docs.conda.io/en/latest/miniconda.html)) environment.

2. Run the following command in a convenient place to install HT-ASE:
```bash
git clone https://github.com/arosen93/htase.git && cd htase && pip install -r requirements.txt && pip install -e .
```

3. Download the example `config` folder provided [here](https://github.com/arosen93/htase/tree/main/htase/setup) and follow the steps in `instructions.md`.

4. Define several environment variables (e.g. in your `~/.bashrc`), as outlined below:
```bash
# HT-ASE requirements
export ASE_VASP_COMMAND="python /path/to/htase/htase/custodian/vasp/run_vasp_custodian.py"
export VASP_CUSTODIAN_SETTINGS="/path/to/config/htase_config/vasp_custodian_settings.yaml"

# ASE requirements
# (details: https://wiki.fysik.dtu.dk/ase/ase/calculators/vasp.html#pseudopotentials)
export VASP_PP_PATH=... # tells ASE where the VASP PAW pseudopotentials are
export ASE_VASP_VDW=... # directory containing vdw_kernel.bindat (optional)

# Jobflow requirements (optional)
# (details: https://materialsproject.github.io/jobflow/jobflow.settings.html)
export JOBFLOW_CONFIG_FILE="/path/to/config/jobflow_config/jobflow.yaml"

# FireWorks requirements (optional)
# (details: https://materialsproject.github.io/fireworks)
export FW_CONFIG_FILE='/path/to/config/fw_config/FW_config.yaml'

```

## Dependencies
### Required
- [ASE](https://gitlab.com/ase/ase)
- [Pymatgen](https://github.com/materialsproject/pymatgen)
- [Custodian](https://github.com/materialsproject/custodian)

### Recommended
- [Jobflow](https://github.com/materialsproject/jobflow)
- [Atomate2](https://github.com/materialsproject/atomate2)
- [FireWorks](https://materialsproject.github.io/fireworks)

## License
HT-ASE is released under a [modified BSD license](https://github.com/arosen93/htase/blob/main/LICENSE.md).


## Representative Examples
### SmartVasp Calculator
For VASP, QuAcc comes pre-packaged with a `SmartVasp` calculator that can run calculations via [Custodian](https://github.com/materialsproject/custodian) for on-the-fly error handling, an optional "co-pilot" mode to ensure your input arguments don't go against what is in the [VASP manual](https://www.vasp.at/wiki/index.php/Main_page), and more.

The example below runs a Custodian-powered VASP relaxation of bulk Cu using the RPBE functional with the remaining settings taken from a [pre-defined](https://github.com/arosen93/quacc/tree/main/quacc/defaults/user_calcs/vasp) set ("preset") of ASE-style calculator input arguments.

```python
from quacc.calculators.vasp import SmartVasp
from ase.build import bulk

atoms = bulk("Cu") # example Atoms object
atoms = SmartVasp(atoms, xc='rpbe', preset="BulkRelaxSet") # set calculator
atoms.get_potential_energy() # run VASP w/ Custodian
```

### Jobflow Integration
The above example can be made compatible with Jobflow simply by defining it in a function with a `@job` wrapper immediately preceeding it. One nuance of Jobflow is that the inputs and outputs must be JSON serializable (so that it can be easily stored in a database), but otherwise it works the same.

```python
from jobflow import job, Job
from quacc.calculators.vasp import SmartVasp
from quacc.schemas.vasp import summarize_run
from quacc.util.json import unclean

#-----Jobflow Function-----
@job
def run_relax(atoms_json):

    # Run VASP
    atoms = SmartVasp(unclean(atoms_json), xc='rpbe', preset="BulkRelaxSet")
    atoms.get_potential_energy()
    
    # Return serialized results
    summary = summarize_run(atoms)
    return summary
```
```python
from ase.build import bulk
from jobflow import Job
from jobflow.managers.local import run_locally
from quacc.util.json import clean

#-----Make and Run a Flow-----
# Constrct an Atoms object
atoms = bulk("Cu")

# Define the job
job1 = run_relax(clean(atoms))

# Run the job locally
run_locally(Job(job1), create_folders=True)
```
### Fireworks Integration
Jobflow provides an easy interface to [Fireworks](https://github.com/materialsproject/fireworks) for high-throughput job management. For additional details on how to convert a Jobflow job or flow to a Fireworks firework or workflow, refer to the [Jobflow documentation](https://materialsproject.github.io/jobflow/jobflow.managers.html#module-jobflow.managers.fireworks). 

### Coupling of Multiple Codes
Since QuAcc is built on top of [ASE](https://wiki.fysik.dtu.dk/ase/index.html) for calculation setup and execution, this means that there is out-of-the-box support for most of your favorite DFT packages. Additionally, through the use of [cclib](https://github.com/cclib/cclib) and [pymatgen](https://pymatgen.org), there's rarely a need to construct your own output parsers. QuAcc will parse most things you can throw at it.

The example below highlights how one could construct a Fireworks workflow to carry out a structure relaxation of O2 using Gaussian and then a static calculation using Q-Chem. The metadata and tabulated calculation results of this workflow would be deposited in your database.
```python
from ase.calculators.gaussian import Gaussian
from ase.calculators.qchem import QChem
from jobflow import job
from quacc.schemas.cclib import summarize_run
from quacc.util.json import unclean

# -----Jobflow Function-----
@job
def run_relax_gaussian(atoms_json):
    # Run Gaussian
    atoms = unclean(atoms_json)
    atoms.calc = Gaussian(method="wB97X-D", basis="def2-TZVP", extra="opt")
    atoms.get_potential_energy()

    # Return serialized results
    summary = summarize_run(atoms)
    return summary


# -----Jobflow Function-----
@job
def run_relax_qchem(atoms_json):

    # Run Q-Chem
    atoms = unclean(atoms_json)
    atoms.calc = QChem(method="wB97M-V", basis="def2-TZVPD", jobtype="SP")
    atoms.get_potential_energy()

    # Return serialized results
    summary = summarize_run(atoms)
    return summary
```
```python
from ase.build import molecule
from fireworks import LaunchPad
from jobflow import Flow
from jobflow.managers.local import flow_to_workflow
from quacc.util.json import clean

# -----Create a Fireworks workflow-----
# Constrct an Atoms object
atoms = molecule("O2")

# Define the flow
job1 = run_relax_gaussian(clean(atoms))
job2 = run_relax_qchem(job1.output["atoms"])
flow = Flow([job1, job2])

# Add a Firework to the launchpad
wf = flow_to_workflow(flow)
lpad = LaunchPad.auto_load()
lpad.add_wf(wf)
```

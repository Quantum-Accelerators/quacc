<div align="center">
  <img src=quacc_logo_wide.svg width="500"><br>
</div>

--------------------------------------
# QuAcc – The Quantum Accelerator
![tests](https://github.com/arosen93/quacc/actions/workflows/tests.yaml/badge.svg)
[![codecov](https://codecov.io/gh/arosen93/quacc/branch/main/graph/badge.svg?token=BCKGTD89H0)](https://codecov.io/gh/arosen93/quacc)
[![CodeFactor](https://www.codefactor.io/repository/github/arosen93/quacc/badge)](https://www.codefactor.io/repository/github/arosen93/quacc)
[![DeepSource](https://deepsource.io/gh/arosen93/quacc.svg/?label=active+issues&token=O0LvluUkUS6qiQnHXc7BUlHn)](https://deepsource.io/gh/arosen93/quacc/?ref=repository-badge)
[![This project supports Python 3.10](https://img.shields.io/badge/Python-3.10+-blue.svg)](https://python.org/downloads)
[![Pypi](https://img.shields.io/pypi/v/quacc)](https://pypi.org/project/quacc)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

QuAcc is a platform for high-throughput, database-driven computational materials science and quantum chemistry. Primarily, QuAcc seeks to enable a seamless interface between the [Atomic Simulation Environment](https://wiki.fysik.dtu.dk/ase/index.html) (ASE) and [Jobflow](https://github.com/materialsproject/jobflow) for rapid workflow development and prototyping while leveraging many of the tools that power the [Materials Project](https://materialsproject.org).

**Disclaimer**: Currently, this package is under active development and should be considered *highly experimental.*

## Documentation
[**Coming Soon**](https://arosen93.github.io/quacc/)


## Examples
### VASP Job
```python
from ase.build import bulk
from jobflow.managers.local import run_locally

from quacc.recipes.vasp.core import RelaxJob as VaspRelaxJob

# Make a bulk Cu structure
atoms = bulk("Cu")

# Make a job consisting of a VASP relaxation using a pre-defined input set.
# By default, VASP will be run using Custodian for on-the-fly error handling.
job = VaspRelaxJob(preset="BulkSet").make(atoms)

# Run the job locally, with all output data stored in a convenient schema
responses = run_locally(job, create_folders=True)
```

### GFN2-xTB + Gaussian + ORCA Workflow
```python
from ase.build import molecule
from fireworks import LaunchPad
from jobflow import Flow
from jobflow.managers.local import run_locally

from quacc.recipes.xtb.core import RelaxJob as XTBRelaxJob
from quacc.recipes.gaussian.core import RelaxJob as GaussianRelaxJob
from quacc.recipes.orca.core import StaticJob as OrcaStaticJob

# Make an H2 molecule
atoms = molecule("H2")

# Make a flow consisting of a GFN2-xTB relaxation followed by a Gaussian relaxation
# and then an ORCA static calculation
job1 = XTBRelaxJob(method="GFN2-xTB").make(atoms)
job2 = GaussianRelaxJob(xc="PBE").make(job1.output["atoms"])
job3 = OrcaStaticJob(xc="wB97M-V").make(job2.output["atoms"])

flow = Flow([job1, job2, job3])

responses = run_locally(flow, create_folders=True)
```

### Database-Friendly Output
Assuming a Jobflow configuration file has been provided, the input and output data will be automagically tabulated and placed in your selected database. No custom parsing required. An example document is shown below:

![docs](docs/src/imgs/schema.gif)

## Inspiration
This package is heavily inspired by [Atomate2](https://github.com/materialsproject/atomate2), which I also recommend checking out.

## License
QuAcc is released under a [modified BSD license](https://github.com/arosen93/quacc/blob/main/LICENSE.md).

# QuAcc (🚧 Under Construction 🚧)
![example workflow](https://github.com/arosen93/quacc/actions/workflows/workflow.yaml/badge.svg)
[![codecov](https://codecov.io/gh/arosen93/quacc/branch/main/graph/badge.svg?token=BCKGTD89H0)](https://codecov.io/gh/arosen93/quacc)
[![CodeFactor](https://www.codefactor.io/repository/github/arosen93/quacc/badge)](https://www.codefactor.io/repository/github/arosen93/quacc)
[![This project supports Python 3.10](https://img.shields.io/badge/Python-3.10+-blue.svg)](https://python.org/downloads)

The Quantum Accelerator (QuAcc) supercharges your code to support high-throughput, database-driven density functional theory (DFT). Primarily, QuAcc seeks to enable a seamless interface between the [Atomic Simulation Environment](https://wiki.fysik.dtu.dk/ase/index.html) (ASE) and [Jobflow](https://github.com/materialsproject/jobflow) for rapid workflow development and prototyping, no matter your favorite DFT package.

This package is heavily inspired by [Atomate2](https://github.com/materialsproject/atomate2), which I also recommend checking out.

## Example: EMT + VASP Flow
```python
from ase.build import bulk
from jobflow import Flow
from jobflow.managers.local import run_locally

from quacc.recipes.emt.core import RelaxMaker as EMTRelaxMaker
from quacc.recipes.vasp.core import RelaxMaker as VaspRelaxMaker

# Make a bulk Cu structure
atoms = bulk("Cu")

# Make a flow consisting of an EMT relaxation followed by a VASP relaxation
job1 = EMTRelaxMaker().make(atoms)
job2 = VaspRelaxMaker(preset="BulkRelaxSet").make(job1.output["atoms"])
flow = Flow([job1, job2])

# Run the flow locally, with all output data stored in a convenient schema
responses = run_locally(flow)
```

## Example: GFN-FF + Gaussian + ORCA Flow with FireWorks
```python
from ase.build import molecule
from fireworks import LaunchPad
from jobflow import Flow
from jobflow.managers.fireworks import flow_to_workflow

from quacc.recipes.xtb.core import RelaxMaker as XTBRelaxMaker
from quacc.recipes.gaussian.core import RelaxMaker as GaussianRelaxMaker
from quacc.recipes.orca.core import StaticMaker as OrcaStaticMaker

# Make an H2 molecule
atoms = molecule("H2")

# Make a flow consisting of a GFN-FF relaxation followed by a Gaussian relaxation
# and then an ORCA static calculation
job1 = XTBRelaxMaker(method="GFN-FF").make(atoms)
job2 = GaussianRelaxMaker().make(job1.output["atoms"])
job3 = OrcaRelaxMaker(xc="wB97M-V").make(job2.output["atoms"])
flow = Flow([job1, job2, job3])

# Convert the flow to a FireWorks workflow and add it to launchpad
# Database-friendly results will be deposited in your JobFlow DB
wf = flow_to_workflow(flow)
lpad = LaunchPad.auto_load()
lpad.add_wf(wf)
```

## Installation
1. Run the following command, ideally in a fresh Python 3.10+ environment: `pip install quacc`. For the most recent development version, instead run `pip install git+https://github.com/arosen93/quacc.git`.

2. Follow the instructions in ASE's [documentation](https://wiki.fysik.dtu.dk/ase/ase/calculators/calculators.html#supported-calculators) for how to set up the ASE calculator(s) you plan to use.

3. Define the following environment variables (e.g. in your `~/.bashrc`) if you wish to use Jobflow and/or Fireworks, in addition to any that you have set in Step 2. Example `.yaml` files are provided [here](https://github.com/arosen93/quacc/tree/main/.config).

```bash
# Jobflow requirements
# (details: https://materialsproject.github.io/jobflow/jobflow.settings.html)
export JOBFLOW_CONFIG_FILE="/path/to/config/jobflow_config/jobflow.yaml"

# FireWorks requirements
# (details: https://materialsproject.github.io/fireworks)
export FW_CONFIG_FILE='/path/to/config/fw_config/FW_config.yaml'

```
## License
QuAcc is released under a [modified BSD license](https://github.com/arosen93/quacc/blob/main/LICENSE.md).

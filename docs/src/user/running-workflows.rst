.. _running_workflows:

=================
Running Workflows
=================

Introduction
============






========
Examples
========
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

```python
from ase.build import molecule
from fireworks import LaunchPad
from jobflow import Flow
from jobflow.managers.fireworks import flow_to_workflow

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

# Instead of running locally, we will run the workflow via Fireworks here.
# The commands below convert the flow to a FireWorks workflow and adds it to
# the launchpad. Database-friendly results will be deposited in your JobFlow DB
wf = flow_to_workflow(flow)
lpad = LaunchPad.auto_load()
lpad.add_wf(wf)
```

============
Contributing
============
To contribute to QuAcc, we recommend installing QuAcc in editable mode so that your changes are reflected every time you run QuAcc. This can be done by using the `-e` flag when running `pip install`. Also, you will need to install the development dependencies via `pip install .[dev]`.

To ensure a consistent formatting style, we use `isort` and `black` to format the code. In the base directory of `quacc`, simply run `isort .` and then `black .` (in that order) to format the code. This will ensure your contributions pass the linting check.

To ensure that any changes you have made do not break unit tests, run `pytest` from the the `quacc/tests` directory.

Any new features that you propose to add to QuAcc must have accompanying unit tests. In certain cases, we have monkeypatched certain functions to change their behavior during testing. For instance, we do not want to run VASP directly during unit tests and have mocked the `atoms.get_potential_energy()` function to always return a value of -1.0 during unit tests. Any mocked functions can be found in the `conftest.py` files of the testing directory.
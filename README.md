<div align="center">
  <img src=docs/src/_static/quacc_logo_wide.svg width="500"><br>
</div>

--------------------------------------

# QuAcc – The Quantum Accelerator

![tests](https://github.com/arosen93/quacc/actions/workflows/tests.yaml/badge.svg)
[![codecov](https://codecov.io/gh/arosen93/quacc/branch/main/graph/badge.svg?token=BCKGTD89H0)](https://codecov.io/gh/arosen93/quacc)
[![DeepSource](https://deepsource.io/gh/arosen93/quacc.svg/?label=active+issues&token=O0LvluUkUS6qiQnHXc7BUlHn)](https://deepsource.io/gh/arosen93/quacc/?ref=repository-badge)
[![This project supports Python 3.8](https://img.shields.io/badge/Python-3.8+-blue.svg)](https://python.org/downloads)
[![Pypi](https://img.shields.io/pypi/v/quacc)](https://pypi.org/project/quacc)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7720998.svg)](https://doi.org/10.5281/zenodo.7720998)

QuAcc is a flexible platform for high-throughput, database-driven computational materials science and quantum chemistry. Built around the [Atomic Simulation Environment](https://wiki.fysik.dtu.dk/ase/index.html) (ASE), Quacc seeks to reduce the barrier for developing and deploying complex, mixed-code workflows across heterogeneous compute environments.

This package was originally inspired by [Atomate2](https://github.com/materialsproject/atomate2), which I also recommend checking out.

**Disclaimer**: Currently, this package is under active development. Questions and bug reports are always welcome!

## Installation

QuAcc can be installed as follows:

```bash
# For the current development version (recommended)
pip install git+https://github.com/arosen93/quacc.git

# For the latest PyPI release
pip install quacc
```

## Documentation

[Click me!](https://arosen93.github.io/quacc/)

## Examples

### VASP Job

```python
from ase.build import bulk

from quacc.recipes.vasp.core import RelaxJob as VaspRelaxJob

# Make a bulk Cu structure
atoms = bulk("Cu")

# Make a job consisting of a VASP relaxation using a pre-defined input set.
# By default, VASP will be run using Custodian for on-the-fly error handling.
output = VaspRelaxJob(atoms, preset="BulkSet")
```

### GFN2-xTB + Gaussian + ORCA Workflow with Covalent

```python
import covalent as ct
from ase.build import molecule

from quacc.recipes.tblite.core import RelaxJob as XTBRelaxJob
from quacc.recipes.gaussian.core import RelaxJob as GaussianRelaxJob
from quacc.recipes.orca.core import StaticJob as OrcaStaticJob

# Make an H2 molecule
atoms = molecule("H2")

# Make a workflow consisting of a GFN2-xTB relaxation followed by a Gaussian relaxation
# and then an ORCA static calculation using Covalent to manage the workflow.
@ct.lattice
def workflow(atoms):
    output1 = ct.electron(XTBRelaxJob)(atoms, method="GFN2-xTB")
    output2 = ct.electron(GaussianRelaxJob)(output1["atoms"], xc="PBE")
    output3 = ct.electron(OrcaStaticJob)(output2["atoms"], xc="wB97M-V")
    return output3

# Run the workflow
output = workflow(atoms)
```

### Database-Friendly Output

Assuming a Jobflow configuration file has been provided, the input and output data will be automagically tabulated and placed in your selected database. No custom parsing required. An example document is shown below:

![docs](docs/src/imgs/schema.gif)

## License

QuAcc is released under a [modified BSD license](https://github.com/arosen93/quacc/blob/main/LICENSE.md).

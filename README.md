<div align="center">
  <img src=docs/src/_static/quacc_logo_wide.svg width="500"><br>
</div>

--------------------------------------

# Quacc – The Quantum Accelerator

![tests](https://github.com/arosen93/quacc/actions/workflows/tests.yaml/badge.svg)
[![codecov](https://codecov.io/gh/arosen93/quacc/branch/main/graph/badge.svg?token=BCKGTD89H0)](https://codecov.io/gh/arosen93/quacc)
[![DeepSource](https://deepsource.io/gh/arosen93/quacc.svg/?label=active+issues&token=O0LvluUkUS6qiQnHXc7BUlHn)](https://deepsource.io/gh/arosen93/quacc/?ref=repository-badge)
[![This project supports Python 3.8](https://img.shields.io/badge/Python-3.8+-blue.svg)](https://python.org/downloads)
[![Pypi](https://img.shields.io/pypi/v/quacc)](https://pypi.org/project/quacc)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7720998.svg)](https://doi.org/10.5281/zenodo.7720998)

Quacc is a flexible platform for high-throughput, database-driven computational materials science and quantum chemistry.

Built around the [Atomic Simulation Environment](https://wiki.fysik.dtu.dk/ase/index.html) (ASE), Quacc seeks to:

1. Reduce the barrier for running complex, mixed-code workflows for molecules and materials across heterogeneous compute environments.

2. Promote rapid workflow development and testing via modern workflow managers, such as [Covalent](https://github.com/AgnostiqHQ/covalent).

3. Enable a seamless interface with much of the software infrastructure powering the [Materials Project](https://materialsproject.org).

In short, Quacc is designed to give the user the building blocks necessary to quickly prototype and construct complex quantum-chemical workflows without sacrificing flexibility.

**Disclaimer**: Currently, this package is under active development. Questions and bug reports are always welcome!

## Installation

Quacc can be installed as follows:

```bash
# For the current development version (recommended)
pip install git+https://github.com/arosen93/quacc.git

# For the latest PyPI release
pip install quacc
```

Then simply run `covalent start` to start the Covalent server and UI.

## Documentation

[Read me!](https://arosen93.github.io/quacc/)

## Minimal Example

```python
import covalent as ct
from ase.build import bulk
from quacc.recipes.emt.slabs import BulkToSlabsFlow

@ct.lattice
def workflow(atoms):
    relaxed_slabs = BulkToSlabsFlow().run(atoms)
    return relaxed_slabs

atoms = bulk("Cu")
dispatch_id = ct.dispatch(workflow)(atoms)
result = ct.get_result(dispatch_id)
```

## License

Quacc is released under a [modified BSD license](https://github.com/arosen93/quacc/blob/main/LICENSE.md).

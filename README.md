<div align="center">
  <img src=docs/src/_static/quacc_logo_wide.svg width="500"><br>
</div>

---

# Quacc – The Quantum Accelerator 🦆

![tests](https://github.com/quantum-accelerators/quacc/actions/workflows/tests.yaml/badge.svg)
[![codecov](https://codecov.io/gh/Quantum-Accelerators/quacc/branch/main/graph/badge.svg?token=OJaOZAH30u)](https://codecov.io/gh/Quantum-Accelerators/quacc)
[![DeepSource](https://app.deepsource.com/gh/Quantum-Accelerators/quacc.svg/?label=active+issues&token=Y1NxOLIuFFEqWdjawIYnJNde)](https://app.deepsource.com/gh/Quantum-Accelerators/quacc/?ref=repository-badge)
![Supported python versions](https://img.shields.io/badge/python-3.9%20%7C%203.10-blue)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7720998.svg)](https://doi.org/10.5281/zenodo.7720998)

Quacc is a flexible platform for high-throughput, database-driven computational materials science and quantum chemistry.

The three main goals of Quacc are to:

1. Reduce the barrier for running complex, mixed-code workflows for molecules and materials across distributed compute environments.

2. Promote rapid workflow development and testing via modern workflow management solutions.

3. Enable a seamless interface between the Atomic Simulation Environment and the software infrastructure powering the [Materials Project](https://materialsproject.org).

## Documentation

<p align="center">
  📖 <a href="https://quantum-accelerators.github.io/quacc/"><b><i>Learn More Here!</i></b></a> 📖
</p>

Includes the following sections:

- 🏎️ [Quick Start](https://quantum-accelerators.github.io/quacc/start/demo.html)
- 🔧 [Installation Guide](https://quantum-accelerators.github.io/quacc/install/install.html)
- 🧠 [User Guide](https://quantum-accelerators.github.io/quacc/user/basics.html)
- 🤝 [Developer Guide](https://quantum-accelerators.github.io/quacc/dev/contributing.html)

## Demonstration ✨

```python
import os
import covalent as ct
from ase.build import bulk

from quacc.recipes.emt.core import relax_job
from quacc.recipes.emt.slabs import BulkToSlabsFlow

# Start the Covalent server
os.system("covalent start")

# Define the workflow and set how to execute
@ct.lattice(executor="local")
def workflow(atoms):

    # Relax a bulk structure
    relaxed_bulk = relax_job(atoms)

    # With the relaxed bulk as input, generate and relax slabs
    relaxed_slabs = BulkToSlabsFlow().run(relaxed_bulk["atoms"])

    return relaxed_slabs

# Make an Atoms object of a bulk Cu structure
atoms = bulk("Cu")

# Dispatch the workflow to the Covalent server
# with the bulk Cu Atoms object as the input
dispatch_id = ct.dispatch(workflow)(atoms)

# Fetch the result from the server
result = ct.get_result(dispatch_id)
print(result)
```

![Demo from Quick Start](docs/src/_static/start/start2.gif)

## Citation

If you use Quacc in your work, please cite it as follows:

- A.S. Rosen, Quacc – The Quantum Accelerator, https://doi.org/10.5281/zenodo.7720998.

## Inspiration

Quacc was originally inspired by [Atomate2](https://github.com/materialsproject/atomate2). Refer to the [FAQ section](https://quantum-accelerators.github.io/quacc/about/faq.html) of the documentation for how these codes differ.

## License ⚖️

Quacc is released under a [BSD 3-Clause license](https://github.com/quantum-accelerators/quacc/blob/main/LICENSE.md).

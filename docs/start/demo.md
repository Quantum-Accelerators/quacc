# Demo

Want to get a feel for what's possible with quacc as fast possible? Here we go!

## Installation

Run the following commands in the terminal, ideally in a fresh Python environment:

```bash
pip install --no-cache-dir https://gitlab.com/ase/ase/-/archive/master/ase-master.zip
pip install quacc[covalent]

quacc set WORKFLOW_ENGINE covalent && covalent start
```

Then open the URL printed in the terminal (typically http://localhost:48008) and run a sample workflow below!

!!! Tip

    Don't want to use Covalent? No problem! Quacc supports a [variety of workflow managers](../user/basics/wflow_overview.md) (or none at all!).

## Demo Workflow

This demo workflow will generate a set of surface slabs from bulk Cu and then run a relaxation and static calculation on each generated slab.

```python
import covalent as ct
from ase.build import bulk
from quacc.recipes.emt.slabs import bulk_to_slabs_flow

# Define the Atoms object
atoms = bulk("Cu")

# Dispatch the workflow
dispatch_id = ct.dispatch(bulk_to_slabs_flow)(atoms)

# Fetch the results
result = ct.get_result(dispatch_id, wait=True)
print(result)
```

![Covalent UI](../images/start/start.gif)

## What Next?

Read through the [User Guide](../user/recipes/recipes_intro.md) to learn more about using quacc! And of course, feel free to explore the calculations you just ran in the Covalent UI.

![Covalent UI](../images/start/ui.jpg)

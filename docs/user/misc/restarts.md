# Restarts

## Parsl

`quacc` supports restarting calculations when using Parsl as the workflow engine. This is achieved through Parsl's memoization system, which caches results based on the function and its input arguments.

### How it Works

When using Parsl, `quacc` automatically registers a special hashing function that creates unique identifiers for ASE `Atoms` objects. This allows Parsl to recognize when the same calculation is being attempted on identical atomic structures. More information about Parsl memoization can be found in the [Parsl documentation](https://parsl.readthedocs.io/en/stable/userguide/checkpoints.html).

!!! Warning

    The hashing excludes the `.info` dictionary and calculator, focusing only on main properties (e.g., positions, numbers, cell, ...). For exact details, see the `ase.Atoms.fromdict` method.

### Usage

To enable restart capabilities in your Parsl workflows, you need to:

1. Configure Parsl's checkpointing behavior.
2. Enable caching for the jobs you wish to checkpoint using `cache=True` in the job decorator.

Here's an example:

```python
import parsl
from ase.build import bulk
from parsl.config import Config
from parsl.utils import get_last_checkpoint

from quacc import job
from quacc.recipes.vasp.core import relax_job
from quacc.wflow_tools.customizers import redecorate

config = Config(...)

parsl.load(config)

# Configure Parsl to checkpoint at the end of each task and every 30 minutes
config.checkpoint_mode = "task_exit"
config.checkpoint_period = "00:30:00"

# Load the last checkpoint, necessary when restarting
config.checkpoint_files = get_last_checkpoint()

# Enable restart capabilities by decorating the job with cache=True
restart_enabled_relax_job = redecorate(relax_job, job(cache=True))

# Create an atoms object and run the calculation
atoms = bulk("Cu", cubic=True)

future = relax_job(atoms)
results = future.result()
```

When restarted, Parsl will compare the hash of the various parameters passed to the job function with the hashes stored in the cache. If a match is found, Parsl will return the cached results instead of recomputing them. Results are only reused if the exact same calculation with identical inputs is attempted.

!!! Warning

    In addition to Parsl's in-built hashing functions for standard types (`list`, `dict`, ...) Quacc provides a custom hashing function for ASE `Atoms` objects. If you need to run jobs with other, complex objects as inputs (e.g., `Phonopy` objects) you may need to register additional hashing functions.

### Technical Details

The hashing implementation can be found in `quacc.atoms.core`:

- `_encode_atoms`: Handles the actual encoding of Atoms objects.
- `get_atoms_id_parsl`: Creates Parsl-compatible byte hashes.
- The hashing function is automatically registered when Parsl is selected as the workflow engine.
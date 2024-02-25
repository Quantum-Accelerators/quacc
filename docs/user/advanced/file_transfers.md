# Transferring Files

## Local File Transfers

### Transfers from a Known File Location

Sometimes, you may want to transfer files between jobs. Many recipes within quacc take an optional keyword argument, `copy_files`, that is a dictionary where the keys are the source folders and the values are the source filenames for the files you wish to have copied to the directory where the calculation is ultimately run.

For instance, if you have the files `CHGCAR` and `WAVECAR` stored in `/my/folder`, then you could ensure that is present in the calculation's working directory:

```python
from pathlib import Path
from ase.build import bulk
from quacc.recipes.vasp.core import relax_job

atoms = bulk("Cu")
relax_job(atoms, copy_files={Path("/my/folder"): ["CHGCAR", "WAVECAR"]})
```

The `copy_files` keyword argument also supports glob patterns, such as `WAVECAR*` to copy all files that start with `WAVECAR` from the source.

### Transfers Between Jobs

Sometimes, however, you may not necessarily know _a priori_ where the source file is. For instance, perhaps you want to copy the file `WAVECAR` from a previous job in your workflow that is stored in a unique directory only determined at runtime. In this scenario, you can still use the `copy_files` keyword argument, but you will need to fetch the prior job's directory.

```python
from pathlib import Path
from ase.build import bulk
from quacc.recipes.vasp.core import relax_job, static_job

atoms = bulk("Cu")
results1 = relax_job(atoms)
static_job(results1["atoms"], copy_files={results1["dir_name"]: ["WAVECAR"]})
```

## Non-Local File Transfers

If using quacc across distributed computing environments, there may not be a guarantee that the file you wish to transfer is even on the same machine as your current job. In this scenario, the best approach is to create a dedicated file staging area on your machine of choice where you can store files you wish to transfer.

For instance, if using a job scheduler like Slurm, you could have a step in your job script that copies files (e.g. via `scp`) from a remote location to a local staging area. Then, you can use the `copy_files` keyword argument to copy the files from the staging area to the calculation's working directory.

Alternatively, several of the supported workflow managers have built-in options to enable file staging.

# Transferring Files

## Local File Transfers

### Transfers from a Known File Location

Sometimes, you may want to transfer files between jobs. Every recipe within quacc takes an optional keyword argument `copy_files` that is a list of absolute filepaths to files you wish to have copied to the directory where the calculation is ultimately run.

For instance, if you have a file `WAVECAR` stored in `/path/to/my/file/stage`, then you could ensure that is present in the calculation's working directory as follows:

```python
from ase.build import bulk
from quacc.recipes.vasp.core import relax_job

atoms = bulk("Cu")
relax_job(atoms, copy_files=["/path/to/my/file/stage/WAVECAR"])
```

### Transfers Between Jobs

Sometimes, however, you may not necessarily know _a priori_ where the source file is. For instance, perhaps you want to copy the file `WAVECAR` from a previous job in your workflow that is stored in a unique directory only determined at runtime. In this scenario, you can still use the `copy_files` keyword argument, but you will need to fetch the prior job's directory. This can be done as follows:

```python
impor tos
from ase.build import bulk
from quacc.recipes.vasp.core import relax_job, static_job

atoms = bulk("Cu")
results1 = relax_job(atoms)
static_job(results1, copy_files=[os.path.join(results1["dir_name"], "WAVECAR")])
```

## Non-Local File Transfers

If using quacc across distributed computing environments, there may not be a guarantee that the file you wish to transfer is even on the same machine as your current job. In this scenario, the best approach is to create a dedicated file staging area on your machine of choice where you can store files you wish to transfer.

For instance, if using a job scheduler like Slurm, you could have a step in your job script that copies files (e.g. via `scp`) from a remote location to a local staging area. Then, you can use the `copy_files` keyword argument to copy the files from the staging area to the calculation's working directory.

Alternatively, several of the supported workflow managers have built-in options to enable file staging. For instance, Covalent has a variety of options for [file transfers](https://docs.covalent.xyz/docs/get-started/quick-start#file-transfers) that can be specified in the `@ct.electron` decorator or pre-instantiated `Electron` object, which can be used to automate file staging before the calculation begins (or after the calculation ends).

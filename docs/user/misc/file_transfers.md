# Transferring Files

## Local File Transfers

### Transfers from a Known File Location

Sometimes, you may want to transfer files between jobs. Many recipes within quacc take an optional keyword argument, `copy_files` for this purpose.

The `copy_files` keyword argument takes two forms. The first form is simply a `SourceDirectory`, which is a `str` or `Path` object that points to the directory containing the files you wish to copy (and decompress) to the runtime directory. For instance, if you have gzipped files from a previous run in `/my/folder` that you want to copy and decompress to your new runtime directory, then you could ensure that this is achieved as follows:

```python
from ase.build import bulk
from quacc.recipes.vasp.core import relax_job

relax_job(atoms, copy_files="/my/folder")
```

The second form is a dictionary where the keys are the source folders and the values are the names of the files you wish to have copied to the directory where the calculation is ultimately run. For instance, if you have the files `CHGCAR.gz` and `WAVECAR.gz` stored in `/my/folder`, then you could ensure that they are present (and decompressed) in the calculation's working directory as follows:

```python
from ase.build import bulk
from quacc.recipes.vasp.core import relax_job

atoms = bulk("Cu")
relax_job(atoms, copy_files={"/my/folder": ["CHGCAR.gz", "WAVECAR.gz"]})
```

You can also copy a specific file from a directory while retaining the directory structure. For instance, if you have the file `/my/folder/WAVECAR.gz` and you want to copy it to the calculation's working directory as `folder/WAVECAR.gz`, you can do so as follows:

```python
relax_job(atoms, copy_files={"/my": ["folder/WAVECAR.gz"]})
```

!!! Tip "Glob Patterns"

    The `copy_files` keyword argument also supports [glob patterns](https://docs.python.org/3/library/pathlib.html#pathlib.Path.glob) for the filenames, such as `WAVECAR*` to copy all files that start with `WAVECAR` from the source.

### Transfers Between Jobs

Sometimes, however, you may not necessarily know _a priori_ where the source file is. For instance, perhaps you want to copy (and decompress) the file `WAVECAR.gz` from a previous job in your workflow that is stored in a unique directory only determined at runtime. In this scenario, you can still use the `copy_files` keyword argument, but you will need to fetch the prior job's directory.

```python
from pathlib import Path
from ase.build import bulk
from quacc.recipes.vasp.core import relax_job, static_job

atoms = bulk("Cu")
results1 = relax_job(atoms)
static_job(results1["atoms"], copy_files={results1["dir_name"]: ["WAVECAR.gz"]})
```

## Non-Local File Transfers

If using quacc across distributed computing environments, there may not be a guarantee that the file you wish to transfer is even on the same machine as your current job. In this scenario, the best approach is to create a dedicated file staging area on your machine of choice where you can store files you wish to transfer.

For instance, if using a job scheduler like Slurm, you could have a step in your job script that copies files (e.g. via `scp`) from a remote location to a local staging area. Then, you can use the `copy_files` keyword argument to copy the files from the staging area to the calculation's working directory.

Alternatively, several of the supported workflow managers have built-in options to enable file staging.

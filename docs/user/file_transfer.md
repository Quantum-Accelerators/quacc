# Transferring Files

Many times, it is necessary to transfer the output files from one job into the working directory of a subequent job. In the context of density functional theory calculations, for instance, this might be passing the converged wavefunction from one calculation to the next so it can be read as an initial guess.

## Within the Job Definition

If you know that a given file is always going to be copied from one job to the next in your workflow, you can specify this when creating your job definition.

Every compute job in quacc uses a utility form [`quacc.util.calc`]() to run calculations, whether it be [`quacc.util.calc.run_calc`](https://quantum-accelerators.github.io/quacc/reference/quacc/util/calc.html#quacc.util.calc.run_calc), [`quacc.util.calc.run_ase_opt`](https://quantum-accelerators.github.io/quacc/reference/quacc/util/calc.html#quacc.util.calc.run_ase_opt), or [`quacc.util.calc.run_ase_vib`](https://quantum-accelerators.github.io/quacc/reference/quacc/util/calc.html#quacc.util.calc.run_ase_vib). Each one of these functions has a keyword argument `copy_files` that takes a list of files to copy from the prior job.

For instance, to copy the `WAVECAR` file from a previous VASP calculation, you could do the following:

```python
from quacc.util.calc import run_calc

run_calc(atoms, copy_files=["WAVECAR"])
```

By default, quacc will know to look for the `WAVECAR` file in the working directory of the prior job since this is actually the current working directory before the new job is run.

## With a Workflow Manager

Sometimes, you may want to copy files between jobs that are not necessarily adjacent in the workflow or in way that is not necessarily part of the job definition itself. For instance, you may want to copy the `WAVECAR` file from a VASP calculation that was run in a previous workflow to a new VASP calculation that is part of a new workflow.

This can be done in a variety of ways depending on your chosen workflow manager.

- If you are using Covalent, refer to the ["How-to Guide"](https://docs.covalent.xyz/docs/user-documentation/how-to/how-to-guide) on transferring files. There are ways of specifying in the `@ct.electron` metadata how to transfer files, both locally and remotely.

- If you are using Parsl, refer to the ["Staging Data Files"](https://parsl.readthedocs.io/en/stable/userguide/data.html#staging-data-files) section of the documentation for how to specify file transfers in the `@python_app` definition.

- If you are using Prefect, there exist a variety of [file transfer options](https://docs.prefect.io/2.10.21/concepts/filesystems/) that can be added as a dedicated task in your workflow.

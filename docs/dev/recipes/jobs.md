# Writing Jobs

## Overview

Individual compute tasks are called jobs. You can think of these as distinct calculations you would submit to the job scheduler (e.g. Slurm). Refer to the [Workflow Syntax](../../user/basics/wflow_syntax.md) page for more information on the syntax of jobs.

## Examples

For an example job using an ASE optimizer, refer to [quacc.recipes.tblite.core.relax_job][].

For a representative job involving a code that uses an internal geometry optimizer, refer to [quacc.recipes.vasp.core.relax_job][].

## Requirements and Recommendations

In quacc, jobs are defined as Python functions that are decorated with `#!Python @job`.

The following are typical guidelines for jobs:

- The function inputs and outputs must be serializable with [pickle](https://docs.python.org/3/library/pickle.html).

- The function generally should take an `Atoms` object as the input, along with any other keyword arguments needed to set calculator parameters and define the logic flow.

- The returned object should typically be a dictionary of results, as defined in `quacc.schemas`. If a dictionary is returned, the keys cannot have "." in them.

- Calculations should be executed using a runner from `quacc.runners`.

- When given the choice between using an ASE optimizer and the electronic structure package's built-in optimizer, you should typically use the latter.

## Schemas

A schema is a dictionary containing tabulated input and output properties from a calculation.

All available schemas are stored in `quacc.schemas`, and type-hints associated with the outputs are stored in `quacc.types`.

All schemas create quacc-compatible dictionaries from a given calculation and also will automatically store results in the user's database, if one is specified in the global settings.

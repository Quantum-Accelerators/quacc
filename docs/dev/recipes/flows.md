# Writing Flows

## Simple Flows

### Overview

A flow is a Python function decorated by `#!Python @flow` that contains a collection of jobs. Refer to the [Workflow Syntax](../../user/basics/wflow_syntax.md) page for more information on the syntax of flows.

The flow itself does not contain computationally intensive tasks. It simply calls other jobs (and/or flows) and defines the workflow logic.

### Production Example

A simple, representative flow can be found in [quacc.recipes.vasp.mp.mp_relax_flow][].

## Dynamic Flows

quacc fully supports complex, dynamic flows where the number of jobs is not known _a priori_. In this case, a common pattern is the use of a subflow, defined with a `#!Python @subflow` decorator. A subflow is just like a flow, except the returned object is a list of job outputs.

!!! Note "Minimal Example"

    Let's do the following:

    1. Add two numbers (e.g. `#!Python 1 + 2`)
    2. Make a list of random length of the output of Step 1 (e.g. `#!Python [3, 3, 3]`)
    3. Add a third number to each element of the list from Step 2 (e.g. `#!Python [3 + 3, 3 + 3, 3 + 3]`)

    In practice, we would want each of the two tasks to be their own compute job.

    ```mermaid
    graph LR
    A[Input] --> B(add) --> C(make list)
    C(make list) --> D(add) --> G[Output]
    C(make list) --> E(add) --> G[Output]
    C(make list) --> F(add) --> G[Output]
    ```

    ```python
    from random import randint
    from quacc import flow, job, subflow

    @job
    def add(a, b):
        return a + b

    @job
    def make_list(val):
        return [val] * randint(2, 5)

    @subflow
    def add_distributed(vals, c):
        outputs = []
        for val in vals:
            output = add(val, c)
            outputs.append(output)
        return outputs

    @flow
    def workflow(a, b, c):
        output1 = add(a, b)
        output2 = make_list(output1)
        output3 = add_distributed(output2, c)
        return output3
    ```

### Production Example

A representative example of a complex flow can be found in [quacc.recipes.emt.slabs.bulk_to_slabs_flow][].

This flow is "dynamic" because the number of slabs is not known until the pre-requisite relaxation step is completed.

!!! Note

    To prevent code duplication, we put common recipes that are code-agnostic in `quacc.recipes.common`. These recipes are not meant to be run directly, but rather called by other recipes.

    For instance, [quacc.recipes.emt.slabs.bulk_to_slabs_flow][] and [quacc.recipes.vasp.slabs.bulk_to_slabs_flow][] both call the same [quacc.recipes.common.slabs.bulk_to_slabs_subflow][].

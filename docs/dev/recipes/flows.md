# Writing Flows

## Simple Flows

### Overview

A flow is a Python function decorated by @flow that contains a collection of jobs. Refer to [](../user/basics/wflow_syntax.md) for more information on the syntax of flows.

The flow itself does not contain computationally intensive tasks. It simply calls other jobs (and/or flows) and defines the workflow logic.

### Examples

A simple, representative flow can be found in [quacc.recipes.vasp.mp.mp_relax_flow][]. More complex flows will be discused after introducing [subflows](subflows.md).

## Dynamic Flows

quacc fully supports complex, dynamic flows where the number of jobs is not known _a priori_.

### Examples

A representative example of a complex flow can be found in [quacc.recipes.emt.slabs.bulk_to_slabs_flow][].

This flow is "dynamic" because the number of slabs is not known until the pre-requisite relaxation step is completed.

!!! Note

    To prevent code duplication, we put common recipes that are code-agnostic in `quacc.recipes.common`. These recipes are not meant to be run directly, but rather called by other recipes.

    For instance, [quacc.recipes.emt.slabs.bulk_to_slabs_flow][] and [quacc.recipes.vasp.slabs.bulk_to_slabs_flow][] both call the same [quacc.recipes.common.slabs.common_bulk_to_slabs_flow][].

### Subflows

A common pattern in complex flows is the use of a subflow, defined with a @subflow decorator. A subflow is just like a flow, except the returned object is a list of job outputs. The need for subflows will become clear in the next section on [futures](futures.md).

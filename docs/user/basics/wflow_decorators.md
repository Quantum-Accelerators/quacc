# Workflow Decorators

## Introduction

Here, we provide code snippets for several decorator-based workflow engines. For a comparison of the pros and cons of each approach, refer to the [Workflow Engines Overview](wflow_overview.md) page. We describe the specifics of how to use each workflow engine in more detail later in the documentation.

## Background

To help enable interoperability between workflow engines, quacc offers a unified set of decorators: [`#!Python @job`](https://quantum-accelerators.github.io/quacc/reference/quacc/wflow_tools/decorators.html#quacc.wflow_tools.decorators.job), [`#!Python @flow`](https://quantum-accelerators.github.io/quacc/reference/quacc/wflow_tools/decorators.html#quacc.wflow_tools.decorators.flow), and [`#!Python @subflow`](https://quantum-accelerators.github.io/quacc/reference/quacc/wflow_tools/decorators.html#quacc.wflow_tools.decorators.subflow).

### `@job`

A `#!Python @job` in quacc is an individual compute task. In the context of an HPC environment, this can typically be thought of as a compute task you would normally submit as a single job to the scheduler.

### `@flow`

A `#!Python @flow` in quacc is a collection of one or more jobs. It defines the logic of the workflow, connecting the output of one job to the inputs of one or more other jobs. As a general guideline, the flow should not contain computations beyond those of the underlying jobs or subflows.

### `@subflow`

A `#!Python @subflow` in quacc is any workflow that returns a list of job outputs and where the number of jobs to be called is not necessarily known until runtime.

=== "Parsl"

    Take a moment to read the Parsl documentation's ["Quick Start"](https://parsl.readthedocs.io/en/stable/quickstart.html) to get a sense of how Parsl works. Namely, you should understand the concept of a [`#!Python python_app`](https://parsl.readthedocs.io/en/stable/1-parsl-introduction.html#Python-Apps) and [`#!Python join_app`](https://parsl.readthedocs.io/en/stable/1-parsl-introduction.html?highlight=join_app#Dynamic-workflows-with-apps-that-generate-other-apps), which describe individual compute tasks and dynamic job tasks, respectively.

    <center>

    | Quacc               | Parsl                  |
    | ------------------- | ---------------------- |
    | `#!Python @job`     | `#!Python @python_app` |
    | `#!Python @flow`    | No effect              |
    | `#!Python @subflow` | `#!Python @join_app`   |

    </center>

=== "Covalent"

    Take a moment to learn about the main [Covalent Concepts](https://docs.covalent.xyz/docs/user-documentation/concepts/concepts-index), namely the [`#!Python @ct.electron`](https://docs.covalent.xyz/docs/user-documentation/concepts/covalent-basics#electron) and [`#!Python @ct.lattice`](https://docs.covalent.xyz/docs/user-documentation/concepts/covalent-basics#lattice) decorators, which describe individual compute tasks and workflows, respectively.

    <center>

    | Quacc               | Covalent                                               |
    | ------------------- | ------------------------------------------------------ |
    | `#!Python @job`     | `#!Python @ct.electron`                                |
    | `#!Python @flow`    | `#!Python @ct.lattice`                                 |
    | `#!Python @subflow` | `#!Python @ct.electron`<br>`#!Python @ct.lattice`</br> |

    </center>

=== "Dask"

    Take a moment to read the Dask Delayed documentation [overview page](https://docs.dask.org/en/stable/delayed.html) to get a sense of how the Dask decorators works and the Dask Distributed [quickstart page](https://distributed.dask.org/en/stable/quickstart.html) to understand how to submit tasks to a Dask cluster. Namely, you should understand the `#!Python @delayed` decorator and how to interface with the `Client`.

    <center>

    | Quacc               | Dask                              |
    | ------------------- | --------------------------------- |
    | `#!Python @job`     | `#!Python @delayed`               |
    | `#!Python @flow`    | No effect                         |
    | `#!Python @subflow` | `#!Python delayed(...).compute()` |

    </center>

=== "Redun"

    Take a moment to read the Redun documentation's [Design Overview page](https://insitro.github.io/redun/design.html) to get a sense of how Redun works. Namely, you should understand the `Task` decorator and how to interface with the `Scheduler`.

    <center>

    | Quacc               | Redun            |
    | ------------------- | ---------------- |
    | `#!Python @job`     | `#!Python @task` |
    | `#!Python @flow`    | `#!Python @task` |
    | `#!Python @subflow` | `#!Python @task` |

    </center>

=== "Jobflow"

    Take a moment to read the Jobflow documentation's [Quick Start](https://materialsproject.github.io/jobflow/tutorials/1-quickstart.html) to get a sense of how Jobflow works. Namely, you should understand the `Job` and `Flow` definitions, which describe individual compute tasks and workflows, respectively.

    <center>

    | Quacc               | Jobflow         |
    | ------------------- | --------------- |
    | `#!Python @job`     | `#!Python @job` |
    | `#!Python @flow`    | N/A             |
    | `#!Python @subflow` | N/A             |

    </center>

    !!! Warning

        Due to the difference in how Jobflow handles workflows compared to other supported workflow engines, any quacc recipes that have been pre-defined with a `#!Python @flow` or `#!Python @subflow` decorator (i.e. have `_flow` in the name) cannot be run directly with Jobflow.

The quacc descriptors are drop-in replacements for the specified workflow engine analogue, which we will use for the remainder of the tutorials. Based on the value for the `WORKFLOW_ENGINE` global variable in your [quacc settings](../settings/settings.md), the appropriate decorator will be automatically selected. If the `WORKFLOW_ENGINE` setting is set to `None` (i.e. `quacc set WORKFLOW_ENGINE None`), the decorators will have no effect on the underlying function.

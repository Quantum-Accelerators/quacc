# Workflow Engines Overview

Everyone's computing needs are different, so we ensured that quacc is interoperable with a variety of modern workflow management tools. There are [300+ workflow management tools](https://workflows.community/systems) out there, so we can't possibly support them all. Instead, we have focused on a select few that adopt a similar decorator-based approach to defining workflows with substantial support for HPC systems.

## Summary

!!! Tip "Recommendations"

    Not sure which to choose? In general, we recommend starting with **Parsl** for most HPC users. For a more feature-rich workflow orchestration platform, we recommend trying **Prefect** or **Jobflow** depending on your needs. Some additional opinions on the matter:

    - Covalent: You want a visual dashboard and are prioritizing the use of distributed compute resources, especially cloud compute.
    - Dask: You already are familiar with the Dask ecosystem and are happy to stick with it.
    - Parsl: You want to run many workflows as fast as possible on one or more job scheduler-based HPC machines.
    - Prefect: You want a visual dashboard with a robust workflow management platform and are familiar with the basic concepts of workflow orchestration.
    - Redun: You are running calculations on AWS.
    - Jobflow: You are familiar with using MongoDB.

=== "Covalent"

    [Covalent](https://github.com/AgnostiqHQ/covalent/) is a user-friendly workflow management solution from the company [Agnostiq](https://agnostiq.ai/).

    Pros:

    - Excellent visual dashboard for job monitoring
    - Easy to use in distributed, heterogeneous compute environments
    - Excellent documentation
    - Automatic and simple database integration
    - The compute nodes do not need to be able to connect to the internet

    Cons:

    - It requires a centralized server to be running continuously in order to manage the workflows unless using [Covalent Cloud](https://www.covalent.xyz/cloud/)
    - Support for job scheduler HPC environments is available but not as robust or performant as other solutions
    - High-security HPC environments may be difficult to access via SSH with the centralized server approach
    - Not as widely used as other workflow management solutions

=== "Dask"

    [Dask](https://www.dask.org/) is a popular parallel computing library for Python. We use [Dask Delayed](https://docs.dask.org/en/stable/delayed.html) for lazy function execution, [Dask Distributed](https://distributed.dask.org/en/stable/) for distributed compute, and (optionally) [Dask-Jobqueue](https://jobqueue.dask.org/en/latest/) for orchestrating the execution on HPC machines.

    Pros:

    - Extremely popular
    - Has native support for running on HPC resources
    - It does not involve a centralized server or network connectivity
    - Supports adaptive scaling of compute resources
    - The dashboard to monitor resource usage is very intuitive

    Cons:

    - If the Dask cluster dies, there is no mechanism to gracefully recover the workflow history
    - Monitoring job progress is more challenging and less detailed than other solutions
    - The documentation, while comprehensive, can be difficult to follow given the various Dask components
    - Calculations cannot be submitted remotely or across disparate compute resources

=== "Parsl"

    [Parsl](https://github.com/Parsl/parsl) is a workflow management solution out of Argonne National Laboratory, the University of Chicago, and the University of Illinois. It is well-adapted for running on virtually any HPC environment with a job scheduler.

    Pros:

    - Extremely configurable and deployable for virtually any HPC environment
    - Quite simple to define the workflows and run them from a Jupyter Notebook
    - Thorough documentation and active user community across academia
    - Well-suited for [pilot jobs](https://en.wikipedia.org/wiki/Pilot_job) and advanced queuing schemes
    - Does not rely on maintaining a centralized server

    Cons:

    - The number of different terms can be slightly overwhelming to those less familiar with HPC
    - Monitoring job progress is more challenging and less detailed than other solutions
    - Debugging failed workflows can be difficult
    - The pilot job model is often a new concept to many HPC users that takes some time to understand

=== "Prefect"

    [Prefect](https://github.com/PrefectHQ/prefect) is a workflow orchestration tool that is popular in the data engineering community. It has an excellent dashboard for monitoring workflows.

    Pros:

    - Quite popular among the data engineering community
    - Excellent web-based dashboard for monitoring workflow progress
    - The free version of Prefect Cloud is reasonably generous
    - Can use advanced queueing schemes to manage workflow execution
    - New features are being added regularly and rapidly

    Cons:

    - For those who are less HPC-savvy, some of the concepts can be quite technical
    - If using Prefect Cloud, the compute nodes must have a network connection
    - The dashboard, while useful for monitoring successes and failures, is not ideal for analyzing results
    - The software is geared more towards data engineering than scientific computing, and that is reflected in the features and documentation

=== "Redun"

    [Redun](https://insitro.github.io/redun/) is a flexible workflow management program developed by [Insitro](https://insitro.com/).

    Pros:

    - Extremely simple syntax for defining workflows
    - Has strong support for task/result caching
    - Useful CLI-based monitoring system
    - Very strong AWS support

    Cons:

    - Currently lacks support for typical HPC job schedulers
    - No user-friendly GUI for job monitoring
    - Does not have a particularly active user community
    - Not updated frequently

=== "Jobflow"

    [Jobflow](https://github.com/materialsproject/jobflow) is developed and maintained by the Materials Project team at Lawrence Berkeley National Laboratory and serves as a seamless interface to [FireWorks](https://github.com/materialsproject/fireworks) or [Jobflow Remote](https://github.com/Matgenix/jobflow-remote) for dispatching and monitoring compute jobs.

    !!! Warning

        Jobflow is not yet compatible with the `#!Python @flow` or `#!Python @subflow` decorators used in many quacc recipes and so should only be used if necessary. See [this issue](https://github.com/Quantum-Accelerators/quacc/issues/1061) to track the progress of this enhancement.

    Pros:

    - Native support for a variety of databases
    - Directly compatible with Atomate2
    - Designed with materials science workflows in mind
    - Actively supported by the Materials Project team

    Cons:

    - Is not fully compatible with all the features of `quacc`
    - Parsing the output of a workflow is not as intuitive as other solutions
    - Defining dynamic workflows with Jobflow's `Response` object can be more complex than other solutions

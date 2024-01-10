# Workflow Engines Overview

Everyone's computing needs are different, so we ensured that quacc is interoperable with a variety of modern workflow management tools. There are [300+ workflow management tools](https://workflows.community/systems) out there, so we can't possibly support them all. Instead, we have focused on a select few that adopt a similar decorator-based approach to defining workflows with substantial support for HPC systems.

## Summary

!!! Tip "Recommendations"

    If you are new to workflow engines, we recommend trying either **Parsl** or **Covalent**. If you have a need for speed and appreciate flexibility, Parsl may be the choice for you. If you are looking for a visual dashboard with an emphasis on distributed compute resources, Covalent may be worth considering.

=== "Parsl"

    [Parsl](https://github.com/Parsl/parsl) is a workflow management solution out of Argonne National Laboratory, the University of Chicago, and the University of Illinois. It is well-adapted for running on virtually any HPC environment with a job scheduler.

    Pros:

    - Extremely configurable and deployable for virtually any HPC environment
    - Quite simple to define the workflows
    - Active community, particularly across academia
    - Well-suited for [pilot jobs](https://en.wikipedia.org/wiki/Pilot_job) and advanced queuing schemes
    - Thorough documentation
    - Does not rely on maintaining a centralized server

    Cons:

    - The number of different terms can be slightly overwhelming to those less familiar with HPC
    - Understanding the various configuration options for your HPC setup can be an initial hurdle
    - Monitoring job progress is more challenging and less detailed than other solutions

=== "Covalent"

    [Covalent](https://github.com/AgnostiqHQ/covalent/) is a workflow management solution from the company [Agnostiq](https://agnostiq.ai/).

    Pros:

    - Best-in-class visual dashboard for job monitoring
    - Easy to use in distributed, heterogeneous compute environments
    - Excellent documentation
    - Automatic and simple database integration
    - The compute nodes do not need to be able to connect to the internet, unlike some of its competitors

    Cons:

    - Not as widely used as other workflow management solutions
    - Only supports standard queuing schemes, limiting throughput
    - It requires a centralized server to be running continuously in order to manage the workflows
    - High-security HPC environments may be difficult to access via SSH with the centralized server approach

=== "Dask"

    [Dask](https://www.dask.org/) is a popular parallel computing library for Python. We use [Dask Delayed](https://docs.dask.org/en/stable/delayed.html) for lazy function execution, [Dask Distributed](https://distributed.dask.org/en/stable/) for distributed compute, and (optionally) [Dask-Jobqueue](https://jobqueue.dask.org/en/latest/) for orchestrating the execution on HPC machines.

    Pros:

    - Extremely popular
    - Has significant support for running on HPC resources
    - It does not involve a centralized server or network connectivity
    - Supports adaptive scaling of compute resources

    Cons:

    - If the Dask cluster dies, there is no mechanism to gracefully recover the workflow history.
    - Monitoring job progress is more challenging and less detailed than other solutions
    - The documentation, while comprehensive, can be difficult to follow given the various Dask components

=== "Redun"

    [Redun](https://insitro.github.io/redun/) is a flexible workflow management program developed by [Insitro](https://insitro.com/).

    Pros:

    - Extremely simple syntax for defining workflows
    - Has strong support for task/result caching
    - Useful console-based monitoring system

    Cons:

    - Currently lacks support for typical HPC job schedulers and platforms other than AWS
    - No user-friendly GUI for job monitoring
    - Less active user community than some other options

=== "Jobflow"

    [Jobflow](https://github.com/materialsproject/jobflow) is developed and maintained by the Materials Project team at Lawrence Berkeley National Laboratory and serves as a seamless interface to [FireWorks](https://github.com/materialsproject/fireworks) or [Jobflow Remote](https://github.com/Matgenix/jobflow-remote) for dispatching and monitoring compute jobs.

    **Jobflow**

    Pros:

    - Native support for a variety of databases
    - Directly compatible with Atomate2
    - Designed with materials science in mind
    - Actively supported by the Materials Project team

    Cons:

    - Is not compatible with the `#!Python @flow` decorator used in some quacc recipes
    - Parsing the output of a workflow is not as intuitive as other solutions
    - Defining dynamic workflows with Jobflow's `Response` object can be more complex than other solutions

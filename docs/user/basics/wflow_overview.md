# Workflow Engines Overview

Everyone's computing needs are different, so we ensured that quacc is interoperable with a variety of modern workflow management tools. There are [300+ workflow management tools](https://workflows.community/systems) out there, so we can't possibly support them all. Instead, we have focused on a select few that adopt a similar decorator-based approach to defining workflows with substantial support for HPC systems.

## Summary

!!! Tip

    If you are new to workflow engines, we recommend trying either **Covalent** or **Parsl**. If you are looking to run high-throughput calculations on an HPC system with a job scheduler, we recommend first trying Parsl. If you are looking for a visual dashboard with an emphasis on distributed compute resources, first try Covalent.

=== "Covalent ⭐"

    [Covalent](https://github.com/AgnostiqHQ/covalent/) is a workflow management solution from the company [Agnostiq](https://agnostiq.ai/).

    Pros:

    - Best-in-class visual dashboard for job monitoring
    - Easy to use in distributed, heterogeneous compute environments
    - Excellent documentation
    - Automatic and simple database integration
    - The compute nodes do not need to be able to connect to the internet

    Cons:

    - Not as widely used as other workflow management solutions
    - It requires a centralized server to be running continuously in order to manage the workflows
    - High-security HPC environments may be difficult to access via SSH with the centralized server approach

=== "Parsl ⭐"

    [Parsl](https://github.com/Parsl/parsl) is a workflow management solution out of Argonne National Laboratory, the University of Chicago, and the University of Illinois. It is well-adapted for running on virtually any HPC environment with a job scheduler.

    Pros:

    - Extremely configurable for virtually any HPC environment
    - Quite simple to define the workflows
    - Active community, particularly across academia
    - Well-suited for [pilot jobs](https://en.wikipedia.org/wiki/Pilot_job) and has near-ideal scaling performance
    - Thorough documentation
    - Does not rely on maintaining a centralized server

    Cons:

    - Understanding the various configuration options for your HPC setup can be an initial hurdle
    - Monitoring job progress is more challenging and less detailed than other solutions
    - The concept of always returning a "future" object can be confusing for new users

=== "Prefect"

    !!! Warning

        Prefect support is currently unavailable until [Issue 1013](https://github.com/Quantum-Accelerators/quacc/issues/1013) is closed.

    [Prefect](https://www.prefect.io/) is a workflow management system that is widely adopted in the data science industry.

    Pros:

    - Very popular in the data science industry with an active community
    - Has a nice dashboard to monitor job progress
    - Supports a variety of job schedulers via `dask-jobqueue`
    - Uses a directed acyclic graph-free model for increased flexibility in workflow definitions

    Cons:

    - Lacks documentation for HPC environments, although it supports them
    - Not practical to use if the compute nodes do not support network connections
    - The dashboard stores data for only a 7 day history by default and does not display the full output of each task
    - Sorting out the details of agents, workers, and queues can be challenging
    - The concept of always returning a "future" object can be confusing for new users

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

    [Jobflow](https://github.com/materialsproject/jobflow) is developed and maintained by the Materials Project team at Lawrence Berkeley National Laboratory and serves as a seamless interface to [FireWorks](https://github.com/materialsproject/fireworks) for dispatching and monitoring compute jobs.

    **Jobflow**

    Pros:

    - Native support for a variety of databases
    - Directly compatible with Atomate2
    - Designed with materials science in mind
    - Actively supported by the Materials Project team

    Cons:

    - Parsing the output of a workflow is not as intuitive as other solutions
    - Defining dynamic workflows with Jobflow's `Response` object can be more complex than other solutions
    - Only supports FireWorks for job dispatching at the moment

    **FireWorks**:

    Pros:

    - Well-suited for a variety of job management approaches
    - Helpful dashboard for monitoring job progress

    Cons:

    - FireWorks documentation can be difficult to navigate without prior experience
    - FireWorks can have a steep learning curve due to its many configuration options
    - The reliance on MongoDB can be challenging for new users and certain HPC environments
    - New features are not planned

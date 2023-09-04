# Workflow Engines Overview

Everyone's computing needs are different, so we ensured that quacc is interoperable with a variety of modern workflow management tools. There are [300+ workflow management tools](https://workflows.community/systems) out there, so we can't possibly support them all. Instead, we have focused on a select few that adopt a similar decorator-based approach to defining workflows with substantial support for HPC systems.

## Choosing a Workflow Engine

### Summary

!!! Tip

    Generally, we highly recommend either **Covalent** or **Parsl** for most users. If you are new to workflow engines or would like a helpful UI to monitor workflows, try Covalent. If you have a need for speed and are savvy with supercomputers, try Parsl.

### Pros and Cons

=== "Covalent"

    [Covalent](https://github.com/AgnostiqHQ/covalent/) is a workflow management solution from the company [Agnostiq](https://agnostiq.ai/).

    Summary: Use Covalent if you are looking for a nice UI and don't mind relying on a long-running server for production calculations.

    Pros:

    - Extremely simple to setup and use, even for complex workflows
    - Best-in-class visual dashboard for job monitoring
    - Easy to use in distributed, heterogeneous compute environments
    - Excellent documentation
    - Automatic and simple database integration
    - The compute nodes do not need to be able to connect to the internet

    Cons:

    - Not as widely used as other workflow management solutions
    - It requires a centralized server to be running continuously in order to manage the workflows
    - High-security HPC environments may be difficult to access via SSH with the centralized server approach
    - Not ideal for large numbers of short-duration jobs on remote HPC machines

=== "Parsl"

    [Parsl](https://github.com/Parsl/parsl) is a workflow management solution out of Argonne National Laboratory, the University of Chicago, and the University of Illinois. It is well-adapted for running on virtually any HPC environment with a job scheduler.

    Summary: Use Parsl if you are looking for the most robust solution for HPC machines and don't mind the lack of a UI.

    Pros:

    - Extremely configurable for virtually any HPC environment
    - Relatively simple to define the workflows
    - Active community, particularly across academia
    - Well-suited for [pilot jobs](https://en.wikipedia.org/wiki/Pilot_job) and has near-ideal scaling performance
    - Thorough documentation
    - Does not rely on maintaining a centralized server of any kind

    Cons:

    - Defining the right configuration options for your desired HPC setup can be an initial hurdle
    - Monitoring job progress is more challenging and less detailed than other solutions
    - The concept of always returning a "future" object can be confusing for new users

=== "Redun"

    [Redun](https://insitro.github.io/redun/) is a flexible workflow management program developed by [Insitro](https://insitro.com/).

    Summary: Use Redun if you are specifically interested in running on AWS or K8s and like a terminal-based monitoring approach.

    Pros:

    - Extremely simple syntax for defining workflows.
    - Has strong support for task/result caching.
    - Supports a variety of compute backends.
    - Useful console-based monitoring system.

    Cons:

    - Currently lacks support for typical HPC job schedulers.
    - No user-friendly UI for job monitoring.
    - Less active user community than some other options.

=== "Jobflow"

    [Jobflow](https://github.com/materialsproject/jobflow) is developed and maintained by the Materials Project team at Lawrence Berkeley National Laboratory and serves as a seamless interface to [FireWorks](https://github.com/materialsproject/fireworks) for dispatching and monitoring compute jobs.

    Summary: Use Jobflow if you want to use compatible software with that used by the Materials Project stack.

    **Jobflow**

    Pros:

    - Simple interface for defining individual jobs and workflows
    - Native support for databases
    - Directly compatible with Atomate2
    - Actively supported by the Materials Project team

    Cons:

    - Parsing the output of a workflow is not as intuitive as other solutions
    - Defining dynamic workflows with Jobflow's `Response` object can be more complex than other solutions
    - Only supports FireWorks for job dispatching at the moment

    **FireWorks**:

    Pros:

    - FireWorks is well-suited for a variety of job management approaches
    - Helpful dashboard for monitoring job progress

    Cons:

    - FireWorks documentation can be difficult to navigate without prior experience
    - FireWorks can have a steep learning curve due to its many configuration options
    - The reliance on MongoDB can be challenging for new users and certain HPC environments
    - New features are not currently planned

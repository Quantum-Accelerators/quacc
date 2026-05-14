# Workflow Engines Overview

Everyone's computing needs are different, so we ensured that quacc is interoperable with a variety of modern workflow management tools. There are [300+ workflow management tools](https://workflows.community/systems) out there, so we can't possibly support them all. Instead, we have focused on a select few that adopt a similar decorator-based approach to defining workflows with substantial support for HPC systems.

## Summary

!!! Tip "Recommendations"

    Not sure which to choose? Generally, we recommend Jobflow for users of academic HPC systems with job schedulers. Otherwise, we recommend Prefect. Some additional opinions on the matter.

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

=== "Ray"

    [Ray](https://www.ray.io/) is a high-performance distributed computing framework originally developed at UC Berkeley. Within quacc, we utilize [Ray Core](https://docs.ray.io/en/latest/ray-core/walkthrough.html) for parallel, distributed task execution, allowing workflows to scale seamlessly from a single laptop to a multi-node cluster.

    [Ray](https://www.ray.io/) is a high-performance distributed computing framework originally developed at UC Berkeley. Within quacc, we utilize [Ray Core](https://docs.ray.io/en/latest/ray-core/walkthrough.html) for parallel, distributed task execution, allowing workflows to scale seamlessly from a single laptop to a multi-node cluster.

    **Pros:**

    - **Widely Adopted & Supported:** Extremely popular with a massive community and active development, backed heavily by the broader ML/AI ecosystem.
    - **Lightweight & Intuitive:** Employs a simple `@ray.remote` decorator to easily parallelize Python functions and manage data objects.
    - **Scalability:** Natively handles both single-node and multi-node distributed execution with remarkably low task-scheduling overhead.
    - **Observability:** Features a robust, built-in web dashboard for real-time monitoring of task progress, resource utilization, and system logs.
    - **HPC Compatibility:** Can be deployed on supercomputers and clusters with provided [HPC deployment guides](https://docs.ray.io/en/latest/cluster/vms/user-guides/community/slurm.html).

    **Cons:**

    - **HPC Integration:** Unlike Parsl or Dask-Jobqueue, Ray does not have native wrappers for traditional HPC job schedulers (e.g., Slurm, PBS). Managing node startup (`ray start`) generally requires manual configuration within batch scripts.
    - **Workflow History:** While the built-in dashboard provides excellent live monitoring, retaining a persistent history of completed workflows requires setting up additional tooling or infrastructure.
    - **Remote Orchestration:** Submitting calculations from a local machine to a remote cluster, or coordinating tasks across disjoint compute resources, requires extra networking setup (such as configuring Ray Client).
    - **Ecosystem Focus:** The broader Ray ecosystem (Tune, Serve, Train) is primarily designed for machine learning workloads, meaning some of the extended features and documentation are less relevant to high-throughput chemistry workflows.

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

    Pros:

    - Native support for a variety of databases
    - Directly compatible with Atomate2
    - Designed with materials science workflows in mind
    - Actively supported by the Materials Project team

    Cons:

    - Requires the use of a database like MongoDB, which may not be widely accessible

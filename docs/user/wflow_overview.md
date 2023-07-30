# Workflow Engines Overview

Everyone's computing needs are different, so we ensured that quacc is interoperable with a variety of modern workflow management tools. There are over [300+ workflow management tools](https://workflows.community/systems) out there, so we can't possibly cover them all. Instead, we have focused on the most popular tools that we have tested and found to be compatible with quacc with minimal additional coding required. The recommended solutions below each use a decorator-based approach to defining workflows, and the basic syntax is very similar across all of them.

## Pros and Cons

=== "Covalent"

    The default suggested workflow management solution is [Covalent](https://github.com/AgnostiqHQ/covalent/) due to its ease of use and helpful dashboard. If you are new to workflow management tools or primarily care about computing across distributed resources, this will likely be the best option for you.

    Pros:

    - Extremely simple to setup and use, even for complex workflows
    - Best-in-class visual dashboard for job monitoring
    - Easy to use in distributed, heterogeneous compute environments
    - Excellent documentation
    - Automatic and simple database integration
    - The compute nodes do not need to be able to connect to the internet

    Cons:

    - Still actively in development
    - Not as widely used as other workflow management solutions
    - For HPC, it currently only supports the Slurm job scheduler out-of-the-box, although custom plugins can be defined
    - High-security HPC environments may be difficult to access via SSH with the centralized server approach

=== "Parsl"

    [Parsl](https://github.com/Parsl/parsl) is a workflow management solution out of Argonne National Laboratory, the University of Chicago, and the University of Illinois. It is well-adapted for running on HPC environments with a job scheduler. Parsl is especially useful for "jobpacking," where instead of each task being an individual (e.g. Slurm) job, a long-running worker will continually pull in and distribute tasks within a single Slurm job until the walltime is reached. This can be useful because it is often easier to move through the queue by requesting a single large job rather than many small jobs.

    Pros:

    - Extremely configurable for virtually any HPC environment
    - Relatively simple to define the workflows
    - Active community, particularly across academia
    - Ideal for jobpacking (e.g. packing many compute tasks into a single Slurm job) and has near-ideal scaling performance
    - Thorough documentation
    - Can be run locally without any server overhead

    Cons:

    - Defining the right configuration options for your desired HPC setup can be an initial hurdle
    - Monitoring job progress is more challenging and less detailed than other solutions
    - Challenging to orchestrate workflows with steps running across heterogeneous resources
    - The concept of always returning a "future" object can be confusing for new users

=== "Jobflow"

    If you are affiliated with the Materials Project team at Lawrence Berkeley National Laboratory, we recommend using [Jobflow](https://github.com/materialsproject/jobflow) to define the workflows coupled with [FireWorks](https://github.com/materialsproject/fireworks) to dispatch them since it widely used there.

    **Jobflow**

    Pros:

    - Simple interface for defining individual jobs and workflows
    - Native support for databases
    - Directly compatible with Atomate2
    - Actively supported by the Materials Project team

    Cons:

    - Parsing the output of a workflow is not as intuitive as other solutions
    - Defining dynamic workflows with Jobflow's `Response` object can be slightly more complex than other solutions
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

# Overview

Everyone's computing needs are different, so we ensured that quacc is interoperable with a variety of modern workflow management tools. There are over [300+](https://workflows.community/systems) workflow management tools out there, so we can't possibly cover them all. Instead, we have focused on the most popular tools that we have tested and found to be compatible with quacc with minimal additional coding required. The recommended solutions below each use a decorator-based approach to defining workflows, and the basic syntax is very similar across all of them.

=== "Covalent"

    The default suggested workflow management solution is [Covalent](https://github.com/AgnostiqHQ/covalent/) due to its ease of use and helpful dashboard. If you are new to workflow management tools or primarily care about computing across distributed resources, this will likely be the best option for you.

    Pros:

    - Extremely simple to setup and use, even for complex workflows
    - Best-in-class visual dashboard for job monitoring
    - Designed with HPC environments in mind but also supports cloud compute
    - Excellent documentation
    - Automatic and simple database integration

    Cons:

    - Still actively in development
    - Not as widely used as other workflow management solutions
    - For HPC, it currently only supports the Slurm job scheduler out-of-the-box, although custom plugins can be defined
    - High-security HPC environments may be difficult to access with the centralized server approach

=== "Parsl"

    [Parsl](https://github.com/Parsl/parsl) is a workflow management solution out of Argonne National Laboratory, the University of Chicago, and the University of Illinois. It is well-adapted for running on HPC environments with a job scheduler. Parsl is especially useful for "jobpacking," where instead of each task being an individual (e.g. Slurm) job, a long-running worker will continually pull in and distribute tasks within a single Slurm job until the walltime is reached. This can be useful because it is often easier to move through the queue by requesting a single large job rather than many small jobs.

    Pros:

    - Extremely configurable for virtually any HPC environment
    - Relatively simple to define the workflows
    - Active community, particularly across academia
    - Ideal for jobpacking and has near-ideal scaling performance
    - Thorough documentation

    Cons:

    - Defining the right configuration options for your desired HPC setup can be a challenge
    - The visual dashboard is not as user friendly as some other solutions
    - Challenging to orchestrate workflows with steps running across heterogeneous resources
    - The concept of always returning a "future" object can be confusing for new users

=== "Jobflow + FireWorks"

    If you are affiliated with the Materials Project team at Lawrence Berkeley National Laboratory, we recommend using [Jobflow](https://github.com/materialsproject/jobflow) to define the workflows coupled with [FireWorks](https://github.com/materialsproject/fireworks) to dispatch them since it widely used there.

    Pros (Jobflow):

    - Simple interface for defining individual jobs and workflows
    - Native support for databases
    - Directly compatible with Atomate2

    Cons (Jobflow):

    - Defining dynamic workflows with Jobflow is slightly more complex than other solutions
    - The strong emphasis on using a database can be a barrier to entry

    Pros (FireWorks):

    - FireWorks is well-suited for a variety of job management approaches
    - Helpful dashboard for monitoring job progress

    Cons (FireWorks):

    - FireWorks documentation can be difficult to navigate without prior experience
    - FireWorks can have a steep learning curve due to its many configuration options
    - The need for a connection to a database can be a challenge in certain HPC environments

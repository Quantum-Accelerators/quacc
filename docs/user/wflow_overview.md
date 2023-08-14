# Workflow Engines Overview

Everyone's computing needs are different, so we ensured that quacc is interoperable with a variety of modern workflow management tools. There are [300+ workflow management tools](https://workflows.community/systems) out there, so we can't possibly support them all. Instead, we have focused on a select few that adopt a similar decorator-based approach to defining workflows with substantial support for HPC systems.

## Choosing a Workflow Engine

### Summary

!!! Tip

    Generally, we recommend either **Covalent** or **Parsl** for most users, depending on your prior experience and goals.

My highly subjective ranking system for various features of the supported workflow engines is summarized below.

|   Workflow Engine   | Ease-of-Use[^1] | GUI[^2] | Flexibility[^3] | HPC Batch[^4] | Cloud[^5] | Pilot Jobs[^6] |
| :-----------------: | :-------------: | :-----: | :-------------: | :-----------: | :-------: | :------------: | --- |
|      Covalent       |       😍        |   😍    |       😀        |      😀       |    😍     |       😢       | 😍  |
| Jobflow + FireWorks |       😐        |   😐    |       😀        |      😍       |    😢     |       😐       | 😐  |
|        Parsl        |       😀        |   ☹️    |       😍        |      😍       |    😐     |       😍       | 😍  |
|       Prefect       |       😐        |   😀    |       😀        |      ☹️       |    😍     |       ☹️       | ☹️  |

    [^1]: How easy is it to learn the syntax, set configuration options, and dispatch a workflow?
    [^2]: How useful is the GUI for monitoring and analyzing workflows?
    [^3]: How flexible is the workflow engine in terms of ways to configure job execution?
    [^4]: How well does the workflow engine operate with HPC batch job schedulers?
    [^5]: How well does the workflow engine operate with cloud compute resources?
    [^6]: How easy is it to use a [pilot job](https://en.wikipedia.org/wiki/Pilot_job) model?

### Pros and Cons

=== "Covalent"

    [Covalent](https://github.com/AgnostiqHQ/covalent/) is a workflow management solution from the company [Agnostiq](https://agnostiq.ai/). It is currently the default option in quacc due to its simplicity for new users and diverse range of execution options.

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
    - It requires a centralized server to be running continuously in order to manage the workflows
    - High-security HPC environments may be difficult to access via SSH with the centralized server approach

=== "Jobflow"

    [Jobflow](https://github.com/materialsproject/jobflow) is developed and maintained by the Materials Project team at Lawrence Berkeley National Laboratory and serves as a seamless interface to [FireWorks](https://github.com/materialsproject/fireworks) for dispatching and monitoring compute jobs.

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

=== "Parsl"

    [Parsl](https://github.com/Parsl/parsl) is a workflow management solution out of Argonne National Laboratory, the University of Chicago, and the University of Illinois. It is well-adapted for running on virtually any HPC environment with a job scheduler.

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
    - Challenging to orchestrate workflows with steps running across heterogeneous resources
    - The concept of always returning a "future" object can be confusing for new users

=== "Prefect"

    [Prefect](https://www.prefect.io/) is a workflow management system that is widely adopted in the data science industry.

    !!! Warning

        Prefect support should be considered experimental at this time.

    Pros:

    - Very popular in the data science industry with an active community
    - Useful dashboard to monitor job progress
    - Supports a variety of job schedulers via `dask-jobqueue`
    - Uses a directed acyclic graph-free model for increased flexibility in workflow definitions

    Cons:

    - Lacks documentation for HPC environments, although it supports them
    - Learning the intricacies of Dask can be challenging for new users
    - Requires the compute nodes to be able to make an outbound network connection
    - The dashboard stores results for only a 7 day history by default and does not store the full output of each task

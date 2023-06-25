# FAQ

**Q: Why is Quacc better than `InsertMyFavoriteCode`?**

A: Better is all in the eye of the beholder. If you are happy with the tools you are using, my goal is not to convince you to switch. That said, if you find that writing high-throughput computational workflows has a huge learning curve or simply takes too much time away from your science, then Quacc may be for you.

**Q: Which workflow manager should I use?**

A: Since there is no single best answer for everyone's needs, we designed Quacc such that it can support a variety of options.

If you are new to workflow management tools, we recommend using [Covalent](https://github.com/AgnostiqHQ/covalent/) due to its ease of use and helpful dashboard.

If you find Covalent isn't for you, we recommend [Parsl](https://github.com/Parsl/parsl) since it is well-documented with many options for distributing work on HPC machines. Parsl is particularly well-suited for a variety of HPC environments common in academia.

If you are affiliated with the Materials Project team at Lawrence Berkeley National Laboratory, we recommend using [Jobflow](https://github.com/materialsproject/jobflow) coupled with [FireWorks](https://github.com/materialsproject/fireworks) since it widely used there.

**Q: Do you have any tips for writing new workflows?**

A: When developing new workflows, I strongly recommend writing an initial draft based on a low-cost ASE calculator like `EMT` or `LJ`, which will allow you to test on your local machine with fast turnaround. Also try to keep things modular. Trying to do too much in one function can make it difficult to debug.

**Q: How is Quacc different than just using the [Atomic Simulation Environment (ASE)](https://wiki.fysik.dtu.dk/ase/)**

A: Quacc is designed to supercharge ASE for the efficient design of high-throughput workflows. Quacc makes it possible to easily design complex ASE workflows that can be run on heterogeneous compute environments with modern workflow managers, such as [Covalent](https://github.com/AgnostiqHQ/covalent) and [Jobflow](https://github.com/materialsproject/jobflow). Quacc also makes it possible to combine ASE with the powerful [Materials Project](https://materialsproject.org/) software stack, such as [Custodian](https://github.com/materialsproject/custodian) for on-the-fly error handling. For the database enthusiasts, Quacc is interoperable with [Maggma](https://github.com/materialsproject/maggma) to store calculation results in a readily queryable format like MongoDB.

**Q: How is Quacc different than [Atomate2](https://github.com/materialsproject/atomate2)?**

A: Atomate2 is a fantastic code that can be used to run libraries of computational materials science workflows. Quacc was originally inspired by Atomate2 and can be thought of as a "sister code" to it. That said, there are some important differences.

In Atomate2, workflows are built around the Materials Project software stack and are fundamentally incompatible with much of the [Atomic Simulation Environment (ASE)](https://wiki.fysik.dtu.dk/ase/). In contrast, Quacc was built around ASE from the ground up while still making it possible to leverage much of the software underlying the Materials Project, such as [Custodian](https://github.com/materialsproject/custodian) for automated error correction. Since Quacc is built around ASE, it is also trivial to add new recipes based on any code that has an existing ASE Calculator even if it is not in Quacc yet.

Atomate2 is built around [Jobflow](https://github.com/materialsproject/jobflow). Quacc is predominantly built around [Covalent](https://github.com/AgnostiqHQ/covalent) first but also supports Jobflow.

Atomate2 is a database-first approach to materials science workflows. Quacc is fully compatible with the use of databases but does not require them, with the goal of reducing the barrier to getting started.

Atomate2 is developed and supported by the [Materials Project](http://materialsproject.org/) team at Lawrence Berkeley National Laboratory. Quacc is developed and maintained by the [Rosen Research Group](https://cbe.princeton.edu/people/andrew-rosen) at Princeton University.

Despite these differences, Quacc and Atomate2 are largely interoperable! Feel free to stitch together workflows from both codes as you see fit using Jobflow!

**Q: I see `@ct.electron` across the code. What if I don't want to use Covalent?**

A: The beauty is that you don't have to use Covalent if you don't want to! By default, the `@ct.electron` wrapper is entirely ignored if you use Quacc without a workflow manager or if you use another workflow manager, such as Jobflow. For this reason, you can use most of the functions within Quacc regardless of your choice of workflow manager.

**Q: I want to contribute a new recipe to Quacc, but it will require additional dependencies. Is that a problem?**

A: Absolutely not a problem! Simply make your dependencies optional in the `pyproject.toml` file and use a {obj}`monty.dev.requires` decorator to ensure your recipes and tests only run if your necessary dependencies are installed. Refer to {obj}`quacc.recipes.tblite.core` and the `[tblite]` extras in the `pyproject.toml` file for an example.

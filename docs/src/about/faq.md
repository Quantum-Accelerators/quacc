# FAQ

**Q: How is Quacc different than just using the [Atomic Simulation Environment (ASE)](https://wiki.fysik.dtu.dk/ase/)**

A: Quacc is designed to supercharge ASE for the efficient design of high-throughput workflows. Quacc makes it possible to easily design complex ASE workflows that can be run on heterogeneous compute environments with modern workflow managers, such as [Covalent](https://github.com/AgnostiqHQ/covalent) and [Jobflow](https://github.com/materialsproject/jobflow). Quacc also makes it possible to combine ASE with the powerful [Materials Project](https://materialsproject.org/) software stack, such as [Custodian](https://github.com/materialsproject/custodian) for on-the-fly error handling. For the database enthusiasts, Quacc is interoperable with [Maggma](https://github.com/materialsproject/maggma) to store calculation results in a readily queryable format like MongoDB.

**Q: How is Quacc different than [Atomate2](https://github.com/materialsproject/atomate2)?**

A: Atomate2 is a fantastic code that is maintained as part of the [Materials Project](https://materialsproject.org/) and can be used to write materials science workflows. Quacc was originally inspired by Atomate2 and can be thought of as a "sister code" to it. That said, there are some important differences.

In Atomate2, workflows are built around the Materials Project software stack and are fundamentally incompatible with much of the [Atomic Simulation Environment (ASE)](https://wiki.fysik.dtu.dk/ase/). In contrast, Quacc was built around ASE from the ground up while still making it possible to leverage much of the software underlying the Materials Project, such as [Custodian](https://github.com/materialsproject/custodian) for automated error correction. Since Quacc is built around ASE, it is also trivial to add new recipes based on any code that has an existing ASE Calculator even if it is not in Quacc yet.

Atomate2 is built around [Jobflow](https://github.com/materialsproject/jobflow). Quacc is predominantly built around [Covalent](https://github.com/AgnostiqHQ/covalent) first but also supports Jobflow.

Atomate2 is a database-first approch to materials science workflows. Quacc supports the use of databases but does not require it, with the goal of reducing the barrier to getting started.

Despite these differences, Quacc and Atomate2 are largely interoperable! Feel free to stitch together workflows from both codes as you see fit!

**Q: Can I contribute to Quacc?**
A: Yes! Please do!

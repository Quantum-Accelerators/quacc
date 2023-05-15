# Frequently Asked Questions

**Q: How is QuAcc different than just using the [Atomic Simulation Environment (ASE)](https://wiki.fysik.dtu.dk/ase/)**

A: QuAcc is designed to supercharge ASE for the efficient design of high-throughput workflows. QuAcc makes it possible to easily design complex ASE workflows that can be run on heterogeneous compute environments with modern workflow managers, such as [Covalent](https://github.com/AgnostiqHQ/covalent). QuAcc also makes it possible to combine ASE with the powerful [Materials Project](https://materialsproject.org/) software stack, such as [Custodian](https://github.com/materialsproject/custodian) for on-the-fly error handling. For the database enthusiasts, QuAcc also is interoperable with [Maggma](https://github.com/materialsproject/maggma) to store calculation results in a readily queryable format.

**Q: How is QuAcc different than [Atomate2](https://github.com/materialsproject/atomate2)?**

A: Atomate2 is a fantastic code that is maintained as part of the [Materials Project](https://materialsproject.org/) and can be used for writing materials science workflows. QuAcc was originally inspired by Atomate2 and can be thought of as a "sister code" to it. There are some important but subtle differences, however.

In Atomate2, workflows are built around the Materials Project software stack and are fundamentally incompatible with much of the [Atomic Simulation Environment (ASE)](https://wiki.fysik.dtu.dk/ase/). In contrast, QuAcc was built around ASE from the ground up while still making it possible to leverage the software underlying the Materials Project, such as [Custodian](https://github.com/materialsproject/custodian) for automated error correction. Since QuAcc is built around ASE, it is also trivial to add new recipes based on any code that has an existing ASE Calculator even if it not in QuAcc yet.

Atomate2 is built around [Jobflow](https://github.com/materialsproject/jobflow). For individual compute tasks, QuAcc is workflow manager-agnostic. For pre-packaged workflows, QuAcc is built around [Covalent](https://github.com/AgnostiqHQ/covalent).

Atomate2 provides prepackaged workflows that can be readily used as-is. In most cases, QuAcc assumes that the user is interested in building their own custom workflows and provides the building blocks to do so.

Despite these differences, QuAcc and Atomate2 are entirely interoperable! Feel free to stitch together workflows from both codes as you see fit!

Q: Can I contribute to QuAcc?
A: Yes! Please do!
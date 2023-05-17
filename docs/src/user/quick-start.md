# The Basics

In Quacc, each code comes with pre-packaged jobs and workflows, which we call recipes for short. This tutorial walks you through how to use these provided recipes to run simple calculations that can be tested out on your local machine.

Once you understand the basics, you should move on to the ["Using Quacc Covalent"](covalent.md) or ["Using Quacc with Jobflow"](jobflow/jobflow.md) guides to learn how to use Quacc with the a workflow manager, which allows you to stich together and run complex Quacc workflows across distributed computing resources.

## Pre-Requisites

If you are not yet familiar with the ASE `Atoms` object, you should read the [ASE tutorial](https://wiki.fysik.dtu.dk/ase/ase/atoms.html), as this is the main object used to represent molecules and structures within Quacc. Additionally, it is worthwhile to be familiar with the basics of what an ASE `Calculator` is, which is described in more detail [here](https://wiki.fysik.dtu.dk/ase/ase/calculators/calculators.html).

If you are coming from the Pymatgen ecosystem, there is a module to convert from Pymatgen `Structure`/`Molecule` objects and the ASE `Atoms` object [here](https://pymatgen.org/pymatgen.io.ase.html).

## A Simple Calculation with EMT

Let's start with a simple example. Here, we will use a cheap calculator based on effective medium theory (EMT) to run a static calculation on a bulk structure of copper, as shown below.

```python
from ase.build import bulk
from quacc.recipes.emt.core import static_job

# Make an Atoms object of a bulk Cu structure
atoms = bulk("Cu")

# Run a static calculation on the Atoms object
result = static_job(atoms)
print(result)
```

Walking through step-by-step, we first defined an `Atoms` object representation of the material we wish to run the calculation on. In this example, we have imported the bulk Cu structure from ASE's predefined library of bulk structures (you could also have read in a `.cif` file using the `ase.io.read` function, for instance).

With the `Atoms` object defined, we then imported a desired recipe and instantiated it. In this case, since we want to use EMT, we can look in `quacc.recipes.emt` to see all the available recipes. We are interested in a static calculation so we imported the `static_job` recipe. We then instantiated and ran the recipe by passing in the `Atoms` object we defined earlier.

## A Simple Calculation with GFN2-xTB

If molecules are more your thing, let's consider another simple example. Here, we want to run a structure relaxation of a water molecule using the semi-empirical quantum mechanics method called [GFN2-xTB](https://doi.org/10.1021/acs.jctc.8b01176). This method is conveniently available in [`tblite`](https://github.com/tblite/tblite), which we will use here. The demonstration below shows how to run this calculation with Quacc.

Note that for this example, you will need to install `tblite`, which can be done with `pip install quacc[tblite]`, as noted in the ["Calculator Setup"](../install/codes.md) of the installation instructions.

```python
from ase.build import molecule

from quacc.recipes.tblite.core import relax_job

# Make an Atoms object of a water molecule
atoms = molecule("H2O")

# Run a relaxation using the GFN2-xTB method
result = relax_job(atoms, method="GFN2-xTB")
print(result)
```

Here, we have imported the `relax_job` recipe from `quacc.recipes.tblite`. Most recipes have several optional keyword arguments that you can specify. In this example, the `method="GFN2-xTB"` keyword indicates that we want to use the GFN2-xTB method, which also happened to be the default value.

## A Pair of Serial Calculations

Now let's return to our bulk Cu example from above and start adding on some complexity. Here, we will use EMT to run a relaxation on the bulk Cu structure and then use the output of this calculation as the input to a static calculation with GFN2-xTB. This example nicely highlights how there are no restrictions in terms of how many codes you can use in a single workflow.

```python
from ase.build import bulk
from quacc.recipes.emt.core import relax_job
from quacc.recipes.tblite.core import static_job

# Make an Atoms object of a bulk Cu structure
atoms = bulk("Cu")

# Run a relaxation calculation with EMT
result1 = relax_job(atoms)

# Run a static calculation with GFN2-xTB
result2 = static_job(result1["atoms"], method="GFN2-xTB")
print(result2)
```

By default, the output of all compute jobs in Quacc is a dictionary summarizing the results of the calculation. It always has a key `"atoms"` that contains a copy of the output `Atoms` object. This can be used to pass structure information between jobs, as shown above.

## A Pair of Parallel Calculations

What if you wanted to run two separate calculations in parallel? Well, you might try to do something like the following, but it will run like any normal Python function wherin `atoms2` will not be run until after the calculation on `atoms1` finished.

```python
from ase.build import bulk, molecule
from quacc.recipes.emt.core import relax_job

atoms1 = bulk("Cu")
atoms2 = molecule("N2")
result1 = relax_job(atoms1)
result2 = relax_job(atoms2)
print([result1, result2])
```

To dispatch calculations in a parallel fashion, you need to use a workflow manager, such as [Covalent](covalent.md). Read on to learn how to define workflows with complex connectivity and how to dispatch them across distributed computing resources.

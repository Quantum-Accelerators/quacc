# The Basics

In Quacc, each code comes with pre-packaged jobs and workflows, which we call recipes for short. This tutorial walks you through how to use these provided recipes to run simple calculations that can be tested out on your local machine.

Once you understand the basics, you should move on to the ["Using Quacc with Covalent"](covalent.md) (recommended) guide to learn how to use Quacc with a workflow manager, which allows you to stich together and run complex Quacc workflows across distributed computing resources. Refer to the ["Alternate Workflow Engines"](advanced/alt_workflows/index.md) section for alternate workflow manager options.

## Pre-Requisites

If you are not yet familiar with the ASE [`Atoms`](https://wiki.fysik.dtu.dk/ase/ase/atoms.html) object, you should read the [ASE tutorial](https://wiki.fysik.dtu.dk/ase/ase/atoms.html), as this is the main object used to represent molecules and structures within Quacc. Additionally, it is worthwhile to be familiar with the basics of an ASE [Calculator](https://wiki.fysik.dtu.dk/ase/ase/calculators/calculators.html).

```{hint}
If you are coming from the Pymatgen ecosystem, you can use the {class}`pymatgen.io.ase.AseAtomsAdaptor` class to convert between Pymatgen {obj}`Structure`/{obj}`Molecule` objects and the ASE `Atoms` object.
```

## A Simple Calculation with EMT

Let's start with a simple example. Here, we will use a cheap calculator based on [effective medium theory (EMT)](<https://doi.org/10.1016/0039-6028(96)00816-3>) to run a structure relaxation on a bulk structure of copper, as shown below.

```python
from ase.build import bulk
from quacc.recipes.emt.core import relax_job

# Make an Atoms object of a bulk Cu structure
atoms = bulk("Cu")

# Run a structure relaxation on the Atoms object
result = relax_job(atoms)
print(result)
```

Walking through step-by-step, we first defined an `Atoms` object representation of the material we wish to run the calculation on. In this example, we have imported the bulk Cu structure from ASE's predefined library of bulk structures.

With the `Atoms` object defined, we then imported a desired recipe and instantiated it. In this case, since we want to use EMT, we can look in {obj}`quacc.recipes.emt` to see all the available recipes. We are interested in doing a structure relaxation, so we imported the {obj}`.emt.core.relax_job` recipe. We then instantiated and ran the recipe by passing in the `Atoms` object we defined earlier.

The recipe output (`result`) is a bit too large to print here; nonetheless, for context, it is a dictionary that has the following primary keys:

`>>> ['atoms', 'atoms_info', 'builder_meta', 'chemsys', 'composition', 'composition_reduced', 'density', 'density_atomic', 'dir_name', 'elements', 'formula_anonymous', 'formula_pretty', 'input_structure', 'name', 'nelements', 'nid', 'nsites', 'parameters', 'parameters_opt', 'results', 'structure', 'symmetry', 'trajectory', 'trajectory_results', 'volume']`

The `"atoms"` key contains a copy of the output `Atoms` object, the `"results"` key contains a dictionary of the results of the calculation, and the `"parameters"` key contains a dictionary of the parameters used in the calculation.

```{hint}
You can make an `Atoms` object from common files like a CIF, XYZ, or POSCAR by using the [`ase.io.read`](https://wiki.fysik.dtu.dk/ase/ase/io/io.html) function. For instance, {obj}`from ase.io import read` followed by `atoms = read("</path/to/file>")`.
```

## A Simple Mixed-Code Workflow

Now let's return to our bulk Cu example from above and start adding on some complexity. Here, we will use EMT to run a relaxation on the bulk Cu structure and then use the output of this calculation as the input to a static calculation with the semi-empirical quantum mechanics method GFN2-xTB as implemented in {obj}`.tblite.core.static_job`. This example highlights how there are no restrictions in terms of how many codes you can use in a single workflow.

```{note}
Some recipes require additional setup. Refer to the [Calculator Setup](../install/codes.md##tblite) section for details. Note that `tblite` is currently available via `pip install` on Linux only.
```

```python
from ase.build import bulk
from quacc.recipes.emt.core import relax_job
from quacc.recipes.tblite.core import static_job

# Make an Atoms object of a bulk Cu structure
atoms = bulk("Cu")

# Run a structure relaxation with EMT
result1 = relax_job(atoms)

# Run a static calculation with GFN2-xTB
result2 = static_job(result1, method="GFN2-xTB")
print(result2)
```

```{hint}
The output of most compute jobs is a dictionary summarizing the results of the calculation. It always has a key `"atoms"` that contains a copy of the output `Atoms` object. The first argument to all recipes accepts either the `Atoms` object directly or a dictionary containing it.
```

What happens if the first job fails, you might ask? Then the code will crash, no results will be stored, and you'd have to start from scratch. That'd be sad, but thankfully this is where using a workflow manager can save the day.

Read on to learn how to define workflows with complex connectivity and how to dispatch them across distributed computing resources.

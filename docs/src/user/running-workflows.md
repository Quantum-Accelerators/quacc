# Running Workflow Recipes

In QuAcc, each code comes with pre-packaged jobs and workflows, which we call recipes for short. This tutorial walks you through how to use these provided recipes to run calculations.

## Pre-Requisites

If you are not yet familiar with the ASE `Atoms` object, you should read the ASE tutorial [here](https://wiki.fysik.dtu.dk/ase/ase/atoms.html), as this is the main object used to represent molecules and structures within QuAcc. If you are coming from the Pymatgen ecosystem, there is a module to interconvert between ASE `Atoms` and Pymatgen `Structure`/`Molecule` objects [here](https://pymatgen.org/pymatgen.io.ase.html). Additionally, it is worthwhile to be familiar with the basics of what an ASE `Calculator` is, which is described in more detail [here](https://wiki.fysik.dtu.dk/ase/ase/calculators/calculators.html).

## Running Jobs

In these examples, we'll primarily use recipes centered around `tblite` because it is a research-grade code that is free to use but can still run on your local machine quite quickly. By default, QuAcc does not install quantum chemistry codes on your behalf. So, as noted in [install.md](../install/install.md), you'll want to run `pip install tblite[ase]` to run these examples.

### Static Job

Let's consider a simple example. We want to run a single-point (i.e. static) calculation on a water molecule using the semi-empirical quantum mechanics method called [GFN2-xTB](https://doi.org/10.1021/acs.jctc.8b01176). This method is conveniently available in several packages, including `tblite`, which we will use here. The demonstration below shows how to run this calculation with QuAcc.

```python
from ase.build import molecule
from jobflow.managers.local import run_locally

from quacc.recipes.tblite.core import StaticJob

# Make an Atoms object of a water molecule
atoms = molecule('H2O')

# Make a job consisting of a GFN2-xTB single-point calculation.
# and attach it to the Atoms object
job = StaticJob(method='GFN2-xTB').make(atoms)

# Run the job locally
response = run_locally(job)

# Fetch the output of the job
output = response[job.uuid][1].output
```

There are several important details to highlight.

First, we need to define an `Atoms` object representation of the molecule (or solid) we wish to run the calculation on. In this example, we have imported the water structure from ASE's predefined library of simple molecules, but you could also read in a `.xyz` file using the `ase.io.read` function.

With the `Atoms` object defined, we then need to import a desired recipe and instantiate it. In this case, since we want to use `tblite`, we can look in `quacc.recipes.tblite` to see all the available recipes. Naturally, we are interested in the `StaticJob` recipe. This recipe takes several optional keyword arguments, one of which is the `method`, which we have set as GFN2-xTB (this is the default and would have been set anyway). Then, we have to attach the `Atoms` object to this job, which we do via the `.make()` command. You can think of the `StaticJob` as defining the keyword arguments for the calculation itself, whereas the `.make()` function defines what the job is acting on.

Next, we run the job locally via the `run_locally()` command, which takes as input the job we previously defined. This command returns a jobflow `Response` object, which is a dictionary that contains the output of all jobs that have been run, each of which is defined a unique key. For instance, the above `response` variable looks like the following (albeit with a different `uuid` key).

```python
{'ad103d4c-f5a4-4fe4-bcc3-27c30befa605': {1: Response(output={'builder_meta': {'emmet_version': '0.51.11', 'pymatgen_version': '2023.3.23', 'pull_request': None, 'database_version': None, 'build_date': '2023-04-11 06:11:04.856427'}, 'charge': 0, 'spin_multiplicity': 1, 'natoms': 3, 'elements': [Element H, Element O], 'nelements': 2, 'nelectrons': 10, 'composition': Comp: H2 O1, 'composition_reduced': Comp: H2 O1, 'formula_alphabetical': 'H2 O1', 'formula_pretty': 'H2O', 'formula_anonymous': 'AB2', 'chemsys': 'H-O', 'symmetry': {'point_group': 'C2v', 'rotation_number': 2.0, 'linear': False, 'tolerance': 0.3, 'eigen_tolerance': 0.01, 'matrix_tolerance': 0.1}, 'molecule': Molecule Summary Site: O (0.0000, 0.0000, 0.1193) Site: H (0.0000, 0.7632, -0.4770) Site: H (0.0000, -0.7632, -0.4770), 'atoms_info': {'_id': '236d061cdafdfe2045157eafa058a2a2'}, 'atoms': Atoms(symbols='OH2', pbc=False, initial_magmoms=...), 'parameters': {'method': 'GFN2-xTB', 'accuracy': 1.0, 'max_iterations': 250, 'electronic_temperature': 300.0, 'cache_api': True, 'verbosity': 1}, 'nid': 'LAPTOP-ROFVQCUO.', 'dir_name': '/mnt/c/Users/asros/AppData/Local/Programs/Microsoft VS Code', 'input_structure': {'builder_meta': {'emmet_version': '0.51.11', 'pymatgen_version': '2023.3.23', 'pull_request': None, 'database_version': None, 'build_date': '2023-04-11 06:11:04.851139'}, 'charge': 0, 'spin_multiplicity': 1, 'natoms': 3, 'elements': [Element H, Element O], 'nelements': 2, 'nelectrons': 10, 'composition': Comp: H2 O1, 'composition_reduced': Comp: H2 O1, 'formula_alphabetical': 'H2 O1', 'formula_pretty': 'H2O', 'formula_anonymous': 'AB2', 'chemsys': 'H-O', 'symmetry': {'point_group': 'C2v', 'rotation_number': 2.0, 'linear': False, 'tolerance': 0.3, 'eigen_tolerance': 0.01, 'matrix_tolerance': 0.1}, 'molecule': Molecule Summary Site: O (0.0000, 0.0000, 0.1193) Site: H (0.0000, 0.7632, -0.4770) Site: H (0.0000, -0.7632, -0.4770), 'atoms_info': {}, 'atoms': Atoms(symbols='OH2', pbc=False, calculator=TBLite(...))}, 'results': {'energy': -137.96777594361674, 'free_energy': -137.96777594361674, 'forces': array([[ 2.45393052e-15, -3.56254749e-15, -7.49527168e-01], [-2.23799030e-15, -1.53506496e-01,  3.74763584e-01], [-2.15940226e-16,  1.53506496e-01,  3.74763584e-01]]), 'charges': array([-0.56105836,  0.28052918,  0.28052918]), 'dipole': array([-7.53260730e-17,  5.14066633e-17, -4.76020974e-01])}, 'name': 'tblite-Static'}, detour=None, addition=None, replace=None, stored_data=None, stop_children=False, stop_jobflow=False)}}
```

As you can see, all the major inputs and outputs are summarized here and can be fetched accordingly. To more easily deal with the output, we can query the `Response` object for the given job we are interested in via `output = response[job.uuid][1].output`. Here, `job.uuid` will automatically retrieve the unique `uuid` for the job we ran, and `[1]` will retrieve the first (and only) job that was run. Finally, we access the output of the job via `.output`. For instance, we can now fetch the single-point energy via `output["results]["energy"]`, which should return -137.96777594361674 eV.

If you completed the "Database Configuration" step in the installation instructions, then you will also see this output schema deposited in your database for future retrieval. Otherwise, it will simply be stored in working memory.

### Relaxation Job

Now let's try a slightly more complicated example. The water molecule we imported previously is not the local minimum energy structure at the GFN2-xTB level of theory, so let's optimize it (i.e. relax it). Like before, we'll look in `quacc.recipes.tblite` to find the appropriate recipe, and we'll notice that there is a suitably named `RelaxJob`. The demonstration is shown below.

```python
from ase.build import molecule
from jobflow.managers.local import run_locally

from quacc.recipes.tblite.core import RelaxJob

# Make an Atoms object of a water molecule
atoms = molecule("H2O")

# Make a job consisting of a GFN2-xTB relaxation calculation.
# and attach it to the Atoms object
job = RelaxJob(method="GFN2-xTB").make(atoms)

# Run the job locally
response = run_locally(job)

# Fetch the output of the job
output = response[job.uuid][1].output
```

As you can see, the QuAcc code is nearly the same! All we had to do was swap out the recipe.

## Running Flows

### Relaxation + Thermochemistry Flow

By now, you have seen how to run individual jobs. Very often, you may want to link jobs together to make a workflow (or a `flow` in Jobflow parlance). For instance, you may want to run a relaxation job followed by a vibrational frequency calculation (including thermodynamic corrections) on the relaxed structure. An example of how to do this is shown below.

```python
from ase.build import molecule
from jobflow import Flow
from jobflow.managers.local import run_locally

from quacc.recipes.tblite.core import RelaxJob, ThermoJob

# Make an Atoms object of a water molecule
atoms = molecule("H2O")

# Define a job consisting of a GFN2-xTB relaxation calculation
# and attach it to the Atoms object
job1 = RelaxJob(method="GFN2-xTB").make(atoms)

# Define a second job consisting of a GFN2-xTB thermochemistry calculation and attach it to the Atoms object from the first job.
# Also pass the computed energy from the first job to the second job.
job2 = ThermoJob(method="GFN2-xTB").make(job1.output["atoms"],energy=job1.output["results"]["energy"])

# Define a flow that stitches together both jobs
flow = Flow([job1, job2])

# Run the flow locally
responses = run_locally(flow)

# Fetch the outputs of each job
output1 = responses[job1.uuid][1].output
output2 = responses[job2.uuid][1].output
```

As you can see above, it is similarly straightforward to stitch together jobs to make a `flow`. We first define both jobs separately. Importantly, we need to make sure that the output `Atoms` object and energy is passed into job2, which we can do by querying the `.output` object. Before running, the `job1.output` will only be a reference, but it will transform into the real output after running.

Then we need to define a `flow` that stitches together both jobs. This is done by passing a list of jobs to the `Flow` constructor. Finally, we run the `flow` locally, and fetch the outputs of each individual job via their unique `uuid` values.

### Combining Multiple Codes

One of the major benefits of QuAcc is that it becomes trivial to combine jobs from multiple different codes. For instance, let's consider a toy example where we want to run a geometry optimization on a bulk Cu structure with effective medium theory (EMT) (using ASE's built-in `EMT` calculator) followed by a static calculation with GFN1-xTB (using the `tblite` calculator). The code is shown below.

```python
from ase.build import bulk
from jobflow import Flow
from jobflow.managers.local import run_locally

from quacc.recipes.emt.core import RelaxJob
from quacc.recipes.tblite.core import StaticJob

# Make an Atoms object of a bulk Cu structure
atoms = bulk("Cu")

# Define a job consisting of an EMT relaxation calculation
# and attach it to the Atoms object
job1 = RelaxJob().make(atoms)

# Define a second job consisting of a GFN1-xTB static calculation
# and attach it to the Atoms object from the first job.
job2 = StaticJob(method="GFN1-xTB").make(job1.output["atoms"])

# Define a flow that stitches together both jobs
flow = Flow([job1, job2])

# Run the flow locally
responses = run_locally(flow)

# Fetch the outputs of each job
output1 = responses[job1.uuid][1].output
output2 = responses[job2.uuid][1].output
```

This should be pretty familiar to you by now! Again, the only real difference is swapping out the recipes you want to use.

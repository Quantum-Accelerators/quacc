# Intro to Jobs

In quacc, each code comes with [pre-packaged jobs and workflows](recipes_list.md), which we collectively call recipes. This tutorial walks you through how to use these provided recipes to run simple calculations that can be tested out on your local machine. We'll start with an introduction to jobs, which are individual calculations that can be stitched together.

## Pre-Requisites

If you are not yet familiar with the ASE [`Atoms`](https://wiki.fysik.dtu.dk/ase/ase/atoms.html) object, you should read the [ASE tutorial](https://wiki.fysik.dtu.dk/ase/ase/atoms.html). Additionally, it is worthwhile to be familiar with the basics of an ASE [Calculator](https://wiki.fysik.dtu.dk/ase/ase/calculators/calculators.html). The ["Open Science with ASE: Core Tutorials"](https://ase-workshop-2023.github.io/tutorial/) is an excellent resource for learning about the basics of ASE.

!!! Note

    Since we are not using a workflow engine for these examples, run the following in the command line:

    ```bash
    quacc set WORKFLOW_ENGINE None
    ```

## Examples

### A Simple Calculation with EMT

```mermaid
graph LR
  A[Input] --> B(EMT Relax) --> C[Output]
```

Let's start with a simple example. Here, we will use a cheap calculator based on effective medium theory (EMT) to run a structure relaxation on a bulk structure of copper. We are interested in doing a structure relaxation, so we will use the [quacc.recipes.emt.core.relax_job][] recipe, as demonstrated below. Most recipes in quacc take an `Atoms` object as the first positional argument and will return a dictionary that summarizes the output.

```python
from ase.build import bulk
from quacc.recipes.emt.core import relax_job

# Make an Atoms object of a bulk Cu structure
atoms = bulk("Cu")

# Run a structure relaxation on the Atoms object
result = relax_job(atoms)
print(result)
```

??? Info "Printed Output"

    ```
    {'atoms': Atoms(symbols='Cu', pbc=True, cell=[[0.0, 1.805, 1.805], [1.805, 0.0, 1.805], [1.805, 1.805, 0.0]]),
    'atoms_info': {'_id': 'd4859270a1a67083343bec0ab783f774'},
    'builder_meta': {'build_date': datetime.datetime(2024, 3, 8, 17, 55, 25, 799465),
                      'emmet_version': '0.78.7',
                      'pymatgen_version': '2024.3.1'},
    'chemsys': 'Cu',
    'composition': Composition('Cu1'),
    'composition_reduced': Composition('Cu1'),
    'converged': True,
    'density': 8.971719800606017,
    'density_atomic': 11.761470249999999,
    'dir_name': '/home/rosen/test/quacc-2024-03-08-17-55-25-774290-64642',
    'elements': [Element Cu],
    'fmax': 0.01,
    'formula_anonymous': 'A',
    'formula_pretty': 'Cu',
    'input_atoms': {'atoms': Atoms(symbols='Cu', pbc=True, cell=[[0.0, 1.805, 1.805], [1.805, 0.0, 1.805], [1.805, 1.805, 0.0]]),
                    'atoms_info': {},
                    'builder_meta': {'build_date': datetime.datetime(2024, 3, 8, 17, 55, 25, 797148),
                                      'emmet_version': '0.78.7',
                                      'pymatgen_version': '2024.3.1'},
                    'chemsys': 'Cu',
                    'composition': Composition('Cu1'),
                    'composition_reduced': Composition('Cu1'),
                    'density': 8.971719800606017,
                    'density_atomic': 11.761470249999999,
                    'elements': [Element Cu],
                    'formula_anonymous': 'A',
                    'formula_pretty': 'Cu',
                    'nelements': 1,
                    'nsites': 1,
                    'symmetry': {'crystal_system': <CrystalSystem.cubic: 'Cubic'>,
                                  'number': 225,
                                  'point_group': 'm-3m',
                                  'symbol': 'Fm-3m',
                                  'symprec': 0.1,
                                  'version': '2.3.1'},
                    'volume': 11.761470249999999},
    'name': 'EMT Relax',
    'nelements': 1,
    'nid': 'rosen.',
    'nsites': 1,
    'nsteps': 0,
    'parameters': {'asap_cutoff': False},
    'parameters_opt': {'alpha': 70.0,
                        'max_steps': 1000,
                        'maxstep': 0.2,
                        'optimizer': 'BFGS',
                        'type': 'optimization'},
    'quacc_version': '0.6.10',
    'results': {'energies': array([-0.00568151]),
                'energy': -0.005681511358581304,
                'forces': array([[-3.14977672e-16,  1.01949224e-15,  3.66995236e-16]]),
                'free_energy': -0.005681511358581304,
                'stress': array([1.36192652e-02, 1.36192652e-02, 1.36192652e-02, 0.00000000e+00,
          7.55159347e-17, 1.47492060e-19])},
    'structure': Structure Summary
    Lattice
        abc : 2.5526554800834367 2.5526554800834367 2.5526554800834367
    angles : 60.00000000000001 60.00000000000001 60.00000000000001
    volume : 11.761470249999999
          A : 0.0 1.805 1.805
          B : 1.805 0.0 1.805
          C : 1.805 1.805 0.0
        pbc : True True True
    PeriodicSite: Cu (0.0, 0.0, 0.0) [0.0, 0.0, 0.0],
    'symmetry': {'crystal_system': <CrystalSystem.cubic: 'Cubic'>,
                  'number': 225,
                  'point_group': 'm-3m',
                  'symbol': 'Fm-3m',
                  'symprec': 0.1,
                  'version': '2.3.1'},
    'trajectory': [Atoms(symbols='Cu', pbc=True, cell=[[0.0, 1.805, 1.805], [1.805, 0.0, 1.805], [1.805, 1.805, 0.0]], calculator=SinglePointCalculator(...))],
    'trajectory_results': [{'energies': array([-0.00568151]),
                            'energy': -0.005681511358581304,
                            'forces': array([[-3.14977672e-16,  1.01949224e-15,  3.66995236e-16]]),
                            'free_energy': -0.005681511358581304,
                            'stress': array([1.36192652e-02, 1.36192652e-02, 1.36192652e-02, 0.00000000e+00,
          7.55159347e-17, 1.47492060e-19])}],
    'volume': 11.761470249999999}
    ```

!!! Note "Loading the JSON Output"

    A serialized version of the print output above is also written to the filesystem with the name `quacc_results.json.gz`. It can be read back in as follows:

    ```python
    from monty.serialization import loadfn

    results = loadfn("quacc_results.json.gz")
    ```

### A Mixed-Code Workflow

```mermaid
graph LR
  A[Input] --> B(EMT Relax) --> C(GFN2-xTB Static) --> D[Output]
```

Now let's add on a bit of complexity. Here, we will use EMT to run a relaxation on the bulk Cu structure and then use the output of this calculation as the input to a static calculation with the semi-empirical quantum mechanics method [GFN2-xTB](https://pubs.acs.org/doi/full/10.1021/acs.jctc.8b01176) as implemented in [quacc.recipes.tblite.core.static_job][].

This example highlights how there are no restrictions in terms of how many codes you can use in a single workflow. It also briefly demonstrates how you can specify optional parameters of a given job, as we will elaborate further in the next example.

!!! Note

    Some codes require additional setup, including `tblite`. Refer to the [Calculator Setup](../../install/codes.md) section for details.

```python
from ase.build import bulk
from quacc.recipes.emt.core import relax_job
from quacc.recipes.tblite.core import static_job

# Make an Atoms object of a bulk Cu structure
atoms = bulk("Cu")

# Run a structure relaxation with EMT
result1 = relax_job(atoms)

# Run a static calculation with GFN2-xTB
result2 = static_job(result1["atoms"], method="GFN2-xTB")  # (1)!
print(result2)
```

1. We have to pass an `Atoms` object to the `static_job` recipe, so we parse it from `result1` output dictionary. We also have specified an optional parameter, `method`. Refer to the [function signature](https://quantum-accelerators.github.io/quacc/reference/quacc/recipes/tblite/core.html#quacc.recipes.tblite.core.static_job) for the `static_job` recipe to see the full list of parameters that can be specified.

??? Info "Printed Output"

    ```
    {'atoms': Atoms(symbols='Cu', pbc=True, cell=[[0.0, 1.805, 1.805], [1.805, 0.0, 1.805], [1.805, 1.805, 0.0]]),
    'atoms_info': {'_id': 'd4859270a1a67083343bec0ab783f774',
                    '_old_ids': ['d4859270a1a67083343bec0ab783f774']},
    'builder_meta': {'build_date': datetime.datetime(2024, 3, 8, 17, 56, 18, 773772),
                      'emmet_version': '0.78.7',
                      'pymatgen_version': '2024.3.1'},
    'chemsys': 'Cu',
    'composition': Composition('Cu1'),
    'composition_reduced': Composition('Cu1'),
    'density': 8.971719800606017,
    'density_atomic': 11.761470249999999,
    'dir_name': '/home/rosen/test/quacc-2024-03-08-17-56-17-983251-36940',
    'elements': [Element Cu],
    'formula_anonymous': 'A',
    'formula_pretty': 'Cu',
    'input_atoms': {'atoms': Atoms(symbols='Cu', pbc=True, cell=[[0.0, 1.805, 1.805], [1.805, 0.0, 1.805], [1.805, 1.805, 0.0]]),
                    'atoms_info': {'_id': 'd4859270a1a67083343bec0ab783f774'},
                    'builder_meta': {'build_date': datetime.datetime(2024, 3, 8, 17, 56, 18, 772151),
                                      'emmet_version': '0.78.7',
                                      'pymatgen_version': '2024.3.1'},
                    'chemsys': 'Cu',
                    'composition': Composition('Cu1'),
                    'composition_reduced': Composition('Cu1'),
                    'density': 8.971719800606017,
                    'density_atomic': 11.761470249999999,
                    'elements': [Element Cu],
                    'formula_anonymous': 'A',
                    'formula_pretty': 'Cu',
                    'nelements': 1,
                    'nsites': 1,
                    'symmetry': {'crystal_system': <CrystalSystem.cubic: 'Cubic'>,
                                  'number': 225,
                                  'point_group': 'm-3m',
                                  'symbol': 'Fm-3m',
                                  'symprec': 0.1,
                                  'version': '2.3.1'},
                    'volume': 11.761470249999999},
    'name': 'TBLite Static',
    'nelements': 1,
    'nid': 'rosen.',
    'nsites': 1,
    'parameters': {'accuracy': 1.0,
                    'cache_api': True,
                    'electronic_temperature': 300.0,
                    'max_iterations': 250,
                    'method': 'GFN2-xTB',
                    'verbosity': 1},
    'quacc_version': '0.6.10',
    'results': {'charges': array([-4.63897809e-11]),
                'dipole': array([-7.99968228e-08,  1.96246646e-07, -5.12600360e-08]),
                'energy': -318.8584605831397,
                'forces': array([[5.44452313e-20, 5.39007790e-19, 2.45003541e-19]]),
                'free_energy': -318.8584605831397,
                'stress': array([14.65181119, 14.65181119, 14.65181119,  1.40704868,  1.01943286,
            0.63181703])},
    'structure': Structure Summary
    Lattice
        abc : 2.5526554800834367 2.5526554800834367 2.5526554800834367
    angles : 60.00000000000001 60.00000000000001 60.00000000000001
    volume : 11.761470249999999
          A : 0.0 1.805 1.805
          B : 1.805 0.0 1.805
          C : 1.805 1.805 0.0
        pbc : True True True
    PeriodicSite: Cu (0.0, 0.0, 0.0) [0.0, 0.0, 0.0],
    'symmetry': {'crystal_system': <CrystalSystem.cubic: 'Cubic'>,
                  'number': 225,
                  'point_group': 'm-3m',
                  'symbol': 'Fm-3m',
                  'symprec': 0.1,
                  'version': '2.3.1'},
    'volume': 11.761470249999999}
    ```

### Modifying Parameters of a Job

As demonstrated in the previous example, each recipe accepts optional keyword arguments to modify the default parameters. To identify the possible parameters a function takes, you'll want to check out the function signature and corresponding documentation, like that for [quacc.recipes.emt.core.relax_job][]. Go ahead; click it! Once you do, you'll see from the docstring that there is one required positional argument for the recipe (the `Atoms` object) and several optional keyword arguments, including those related to the underlying ASE calculator (`**calc_kwargs`), the optimization settings, and more. You can modify these parameters like in the example below:

```python
from ase.build import bulk
from ase.optimize import LBFGS
from quacc.recipes.emt.core import relax_job

# Make an Atoms object of a bulk Cu structure
atoms = bulk("Cu")

# Run a structure relaxation on the Atoms object
result = relax_job(
    atoms,
    relax_cell=True,
    opt_params={"fmax": 1e-3, "optimizer": LBFGS},
    asap_cutoff=True,  # (1)!
)
print(result)
```

1. All recipes in quacc allow you to pass in a custom set of keyword arguments corresponding to the ASE calculator.

??? Info "Printed Output"

    ```text
    {'atoms': Atoms(symbols='Cu', pbc=True, cell=[[1.6657488637835028e-19, 1.7963940854299556, 1.7963940854299547], [1.7963940854299552, 2.3931154361445102e-17, 1.7963940854299547], [1.7963940854299552, 1.7963940854299556, -8.132837144101245e-16]]),
    'atoms_info': {'_id': '6915f7a4e111b833dcca453da6d0d694'},
    'builder_meta': {'build_date': datetime.datetime(2024, 3, 8, 17, 54, 26, 361990),
                      'emmet_version': '0.78.7',
                      'pymatgen_version': '2024.3.1'},
    'chemsys': 'Cu',
    'composition': Composition('Cu1'),
    'composition_reduced': Composition('Cu1'),
    'converged': True,
    'density': 9.10127989823884,
    'density_atomic': 11.594041355280432,
    'dir_name': '/home/rosen/test/quacc-2024-03-08-17-54-26-278865-61885',
    'elements': [Element Cu],
    'fmax': 0.001,
    'formula_anonymous': 'A',
    'formula_pretty': 'Cu',
    'input_atoms': {'atoms': Atoms(symbols='Cu', pbc=True, cell=[[0.0, 1.805, 1.805], [1.805, 0.0, 1.805], [1.805, 1.805, 0.0]]),
                    'atoms_info': {},
                    'builder_meta': {'build_date': datetime.datetime(2024, 3, 8, 17, 54, 26, 350899),
                                      'emmet_version': '0.78.7',
                                      'pymatgen_version': '2024.3.1'},
                    'chemsys': 'Cu',
                    'composition': Composition('Cu1'),
                    'composition_reduced': Composition('Cu1'),
                    'density': 8.971719800606017,
                    'density_atomic': 11.761470249999999,
                    'elements': [Element Cu],
                    'formula_anonymous': 'A',
                    'formula_pretty': 'Cu',
                    'nelements': 1,
                    'nsites': 1,
                    'symmetry': {'crystal_system': <CrystalSystem.cubic: 'Cubic'>,
                                  'number': 225,
                                  'point_group': 'm-3m',
                                  'symbol': 'Fm-3m',
                                  'symprec': 0.1,
                                  'version': '2.3.1'},
                    'volume': 11.761470249999999},
    'name': 'EMT Relax',
    'nelements': 1,
    'nid': 'rosen.',
    'nsites': 1,
    'nsteps': 3,
    'parameters': {'asap_cutoff': True},
    'parameters_opt': {'max_steps': 1000,
                        'maxstep': 0.2,
                        'optimizer': 'LBFGS',
                        'type': 'optimization'},
    'quacc_version': '0.6.10',
    'results': {'energies': array([-0.00159299]),
                'energy': -0.0015929904219369462,
                'forces': array([[-1.64798730e-16,  3.27862737e-16,  3.78169718e-15]]),
                'free_energy': -0.0015929904219369462,
                'stress': array([ 1.12813866e-06,  1.12813866e-06,  1.12813866e-06, -1.53212912e-16,
          -2.29819368e-16, -4.92256329e-17])},
    'structure': Structure Summary
    Lattice
        abc : 2.540484878981855 2.5404848789818546 2.540484878981855
    angles : 60.00000000000002 60.00000000000001 60.00000000000002
    volume : 11.594041355280432
          A : 1.6657488637835028e-19 1.7963940854299556 1.7963940854299547
          B : 1.7963940854299552 2.3931154361445102e-17 1.7963940854299547
          C : 1.7963940854299552 1.7963940854299556 -8.132837144101245e-16
        pbc : True True True
    PeriodicSite: Cu (-6.919e-18, 9.514e-17, 1.138e-16) [6.008e-17, 3.264e-18, -7.116e-18],
    'symmetry': {'crystal_system': <CrystalSystem.cubic: 'Cubic'>,
                  'number': 225,
                  'point_group': 'm-3m',
                  'symbol': 'Fm-3m',
                  'symprec': 0.1,
                  'version': '2.3.1'},
    'trajectory': [Atoms(symbols='Cu', pbc=True, cell=[[0.0, 1.805, 1.805], [1.805, 0.0, 1.805], [1.805, 1.805, 0.0]], calculator=SinglePointCalculator(...)),
                    Atoms(symbols='Cu', pbc=True, cell=[[0.0, 1.8014556448273662, 1.8014556448273664], [1.8014556448273662, 0.0, 1.8014556448273664], [1.8014556448273662, 1.8014556448273662, 0.0]], calculator=SinglePointCalculator(...)),
                    Atoms(symbols='Cu', pbc=True, cell=[[5.771123678574472e-17, 1.7963098069799175, 1.7963098069799164], [1.796309806979917, 1.7488253571435684e-18, 1.7963098069799164], [1.796309806979917, 1.7963098069799175, -1.031316315430271e-15]], calculator=SinglePointCalculator(...)),
                    Atoms(symbols='Cu', pbc=True, cell=[[1.6657488637835028e-19, 1.7963940854299556, 1.7963940854299547], [1.7963940854299552, 2.3931154361445102e-17, 1.7963940854299547], [1.7963940854299552, 1.7963940854299556, -8.132837144101245e-16]], calculator=SinglePointCalculator(...))],
    'trajectory_results': [{'energies': array([-0.00060116]),
                            'energy': -0.0006011628239370737,
                            'forces': array([[1.18655086e-15, 1.15706056e-15, 3.16066617e-15]]),
                            'free_energy': -0.0006011628239370737,
                            'stress': array([0.01169831, 0.01169831, 0.01169831, 0.        , 0.        ,
          0.        ])},
                            {'energies': array([-0.00124764]),
                            'energy': -0.0012476410887956035,
                            'forces': array([[-9.64506253e-16,  1.88737914e-15, -8.84708973e-16]]),
                            'free_energy': -0.0012476410887956035,
                            'stress': array([ 6.97489804e-03,  6.97489804e-03,  6.97489804e-03,  0.00000000e+00,
          -7.59625440e-17, -2.37382950e-18])},
                            {'energies': array([-0.0015929]),
                            'energy': -0.0015928956483097778,
                            'forces': array([[6.81746326e-16, 9.74914593e-16, 7.14012183e-15]]),
                            'free_energy': -0.0015928956483097778,
                            'stress': array([-1.17305499e-04, -1.17305499e-04, -1.17305499e-04, -1.53234478e-16,
          -7.66172391e-17,  7.73654544e-17])},
                            {'energies': array([-0.00159299]),
                            'energy': -0.0015929904219369462,
                            'forces': array([[-1.64798730e-16,  3.27862737e-16,  3.78169718e-15]]),
                            'free_energy': -0.0015929904219369462,
                            'stress': array([ 1.12813866e-06,  1.12813866e-06,  1.12813866e-06, -1.53212912e-16,
          -2.29819368e-16, -4.92256329e-17])}],
    'volume': 11.594041355280432}
    ```

# Understanding the Code

## Pre-Requisites

The concepts in this section are based around [Jobflow](https://materialsproject.github.io/jobflow/), namely the concepts of a `job`, `flow`, `Maker`, and `Response`. It may be helpful to review the Jobflow tutorials if you are looking for more detailed information about any of these asepcts used in the examples below.

## Worked Examples

Before writing your own workflows, let's first break down a few examples already provided by QuAcc so that you understand the nuts and bolts.

### Example 1

For this example, I'll walk you through the `quacc.recipes.tblite.core.StaticJob` job step-by-step.

After the somewhat self-explanatory import statements, the first real lines are as follows:

```python
@dataclass
class StaticJob(Maker):
```

Here, we have defined a job with the name `StaticJob`. Strictly speaking, this what is called a jobflow `Maker` class, which is a class that can be used to generate a `job`. You always start by inheriting from the `Maker` class in Jobflow. The `@dataclass` decorator is allows us to define a class with attributes that are automatically initialized when the class is instantiated and should also always be included. So, really the only thing of note here is the name of the class, which is `StaticJob`.

Building it up a bit further, we have

```python
@dataclass
class StaticJob(Maker):
    """
    Class to carry out a single-point calculation.

    Parameters
    ----------
    name
        Name of the job.
    method
        GFN1-xTB, GFN2-xTB, and IPEA1-xTB.
    tblite_kwargs
        Dictionary of custom kwargs for the tblite calculator.
    """

    name: str = "tblite-Static"
    method: str = "GFN2-xTB"
    tblite_kwargs: Dict[str, Any] = field(default_factory=dict)
```

Below the name of the class, we have some docstrings to highlight the inputs that will define the settings of the job to be run. The first argument must always be `name` and is the name of the job. This is the only argument that needs to be included when writing a new `Maker`, but you can include as many additional arguments as you like. Here, we also include an optional argument `method` that will pass to the `tblite` calculator what level of theory to run. Finally, we have an optional argument `tblite_kwargs` that allows the user to pass a dictionary of any additional keyword arguments to the `tblite` calculator. It is important to give the user flexibility, so it is always a good idea to include a dictionary of arguments the user can pass on to the calculator if they so choose.

In case you are not familiar with the syntax here, `name: str = "tblite-Static"` means that the `name` argument takes an argument in the form of a string and by default, this is set to `"tblite-Static"`. The line `tblite_kwargs: Dict[str, Any] = field(default_factory=dict)` means that the `tblite_kwargs` optional argument takes a dictionary of any type of value and by default, this is set to an empty dictionary. Of course, if you are coding this up only for internal use, you can hard-code the default values for the various arguments if you'd like, but this is good programming practice.

Okay, let's continue and expand on our codeblock:

```python
@dataclass
class StaticJob(Maker):
    """
    Class to carry out a single-point calculation.

    Parameters
    ----------
    name
        Name of the job.
    method
        GFN1-xTB, GFN2-xTB, and IPEA1-xTB.
    tblite_kwargs
        Dictionary of custom kwargs for the tblite calculator.
    """

    name: str = "tblite-Static"
    method: str = "GFN2-xTB"
    tblite_kwargs: Dict[str, Any] = field(default_factory=dict)

    @job
    @requires(
        TBLite,
        "tblite must be installed. Try pip install tblite[ase]",
    )
    def make(self, atoms: Atoms) -> Dict[str, Any]:
```

It's time to talk about the `make` function. This is the function that will define the job itself using the aforementioned arguments. The `make` function must always take at least the following arguments: `self` and `atoms`. The former allows the job to inherit the arguments we defined above (i.e. `name`, `method`, `tblite_kwargs`). The latter is necessary so the job knows what `Atoms` object to act on. Again, if you aren't familiar with the syntax here, `atoms: Atoms` means that it is a required positional argument with type `Atoms` and has no default value. Additionally, `-> Dict[str, Any]` is a type hint that tells the user what the output will be. In this case, it is a dictionary of strings and any type of value. Type hints don't change the function of the code; it is just helpful for the user to know what to expect.

The `@job` decorator is equally important! It tells Jobflow that the `make` function is a job and that it should be run. All `job` objects need to have a `@job` wrapper preceeding it. The `@requires` decorator, on the other hand, is never really needed. Here, it is included as a guide for the user. It will throw an error if `tblite` isn't installed, which can be helpful.

Alright, let's wrap this up:

```python
@dataclass
class StaticJob(Maker):
    """
    Class to carry out a single-point calculation.

    Parameters
    ----------
    name
        Name of the job.
    method
        GFN1-xTB, GFN2-xTB, and IPEA1-xTB.
    tblite_kwargs
        Dictionary of custom kwargs for the tblite calculator.
    """

    name: str = "tblite-Static"
    method: str = "GFN2-xTB"
    tblite_kwargs: Dict[str, Any] = field(default_factory=dict)

    @job
    @requires(
        TBLite,
        "tblite must be installed. Try pip install tblite[ase]",
    )
    def make(self, atoms: Atoms) -> Dict[str, Any]:
        """
        Make the run.

        Parameters
        ----------
        atoms
            .Atoms object

        Returns
        -------
        Dict
            Summary of the run.
        """
        atoms.calc = TBLite(method=self.method, **self.tblite_kwargs)
        new_atoms = run_calc(atoms)
        summary = summarize_run(
            new_atoms, input_atoms=atoms, additional_fields={"name": self.name}
        )

        return summary
```

That's the full function! Below `make`, we have some docstrings again highlighting both the inputs and the outputs. Then we do four things: 

1. We instantiate the `TBLite` [ASE calculator](https://tblite.readthedocs.io/en/latest/users/ase.html) with the `method` and any `tblite_kwargs` arguments we defined above. The `**` before th `self.tblite_kwargs` is necessary to unpack the dictionary into individual keyword arguments, if there are any. We also attached this calculator to the `Atoms` object supplied by the user.

2. We then run the calculation using the `run_calc` function supplied with QuAcc. This function is a wrapper around ASE's `.get_potential_energy()` function that will effectively launch a calculation given an `Atoms` object with an attached calculator. Crucially, QuAcc does *not* modify `Atoms` objects in place whenever possible, so we need to save the new `Atoms` object with a new variable name (here, `new_atoms` instead of `atoms`).

3. We create a summary of the run inputs and outputs, which is a dictionary that will be output and stored. This is done using the `summarize_run` function supplied with QuAcc, which takes the output `Atoms`, input `Atoms`, and any other additional fields you'd like to specify. It is always recommended to store the `name` of the job as a field in the summary.

4. We return the summary dictionary to the user.

That's it! Again, from the user perspective you would run this as `StaticJob().make(atoms)` once imported.

### Example 2

Let's continue with the `tblite` examples!

This time, let's take a deep dive with the `quacc.recipes.tblite.core.RelaxJob`.

I'll just give you the whole codeblock to start, but we'll still go through it together.

```python
@dataclass
class RelaxJob(Maker):
    """
    Class to relax a structure.

    Parameters
    ----------
    name
        Name of the job.
    method
        GFN0-xTB, GFN1-xTB, GFN2-xTB.
    fmax
        Tolerance for the force convergence (in eV/A).
    max_steps
        Maximum number of steps to take.
    optimizer
        .Optimizer class to use for the relaxation.
    tblite_kwargs
        Dictionary of custom kwargs for the tblite calculator.
    opt_kwargs
        Dictionary of kwargs for the optimizer.
    """

    name: str = "tblite-Relax"
    method: str = "GFN2-xTB"
    fmax: float = 0.01
    max_steps: int = 1000
    optimizer: str = "FIRE"
    tblite_kwargs: Dict[str, Any] = field(default_factory=dict)
    opt_kwargs: Dict[str, Any] = field(default_factory=dict)

    @job
    @requires(
        TBLite,
        "tblite must be installed. Try pip install tblite[ase]",
    )
    def make(self, atoms: Atoms) -> Dict[str, Any]:
        """
        Make the run.

        Parameters
        ----------
        atoms
            .Atoms object

        Returns
        -------
        Dict
            Summary of the run.
        """
        atoms.calc = TBLite(method=self.method, **self.tblite_kwargs)
        traj = run_ase_opt(
            atoms,
            fmax=self.fmax,
            max_steps=self.max_steps,
            optimizer=self.optimizer,
            opt_kwargs=self.opt_kwargs,
        )
        summary = summarize_opt_run(
            traj, atoms.calc.parameters, additional_fields={"name": self.name}
        )

        return summary
```

So, what's different about this one? Well, the main thing is that some ASE calculators (including `tblite`) require you to use an ASE-based optimizer to relax the structure (see [here](https://wiki.fysik.dtu.dk/ase/ase/optimize.html) for more details). Just like how `run_calc` was a wrapper around ASE's `.get_potential_energy()` function to run a calculation, the `run_ase_opt` function is a wrapper around ASE's `Optimize` class. Otherwise, everything is largely the same except for some different arguments and a different summary function.

### Example 3

Finally, let's look at a more typical relaxation job where you can run the relaxation using the code itself rather than ASE. For this example, we'll look at the `quacc.rescipes.gaussian.core.RelaxJob` class, which is based around the [ASE Gaussian calculator](https://wiki.fysik.dtu.dk/ase/ase/calculators/gaussian.html).

```python

@dataclass
class RelaxJob(Maker):
    """
    Class to carry out a geometry optimization.

    Parameters
    ----------
    name
        Name of the job.
    xc
        Exchange-correlation functional
    basis
        Basis set
    freq
        If a requency calculation should be carried out.
    swaps
        Dictionary of custom kwargs for the calculator.
    """

    name: str = "Gaussian-Relax"
    xc: str = "wb97x-d"
    basis: str = "def2-tzvp"
    freq: bool = False
    swaps: Dict[str, Any] = field(default_factory=dict)

    @job
    def make(
        self, atoms: Atoms, charge: int = None, mult: int = None
    ) -> Dict[str, Any]:
        """
        Make the run.

        Parameters
        ----------
        atoms
            .Atoms object
        charge
            Charge of the system. If None, this is determined from the sum of
            atoms.get_initial_charges().
        mult
            Multiplicity of the system. If None, this is determined from 1+ the sum
            of atoms.get_initial_magnetic_moments().

        Returns
        -------
        Dict
            Summary of the run.
        """
        defaults = {
            "mem": "16GB",
            "chk": "Gaussian.chk",
            "nprocshared": multiprocessing.cpu_count(),
            "xc": self.xc,
            "basis": self.basis,
            "charge": charge,
            "mult": mult,
            "opt": "",
            "scf": ["maxcycle=250", "xqc"],
            "integral": "ultrafine",
            "nosymmetry": "",
            "freq": "" if self.freq else None,
        }
        flags = merge_dicts(defaults, self.swaps, remove_none=True)

        atoms.calc = Gaussian(**flags)
        atoms = run_calc(atoms, geom_file=GEOM_FILE)
        summary = summarize_run(atoms, LOG_FILE, additional_fields={"name": self.name})

        return summary
```

Again, this should largely look similar. We start, like always, by instantiating a class that inherits a jobflow `Maker`. We then define the parameters that we want to be able to specify when we run the job. There are few more parameters here than in the `tblite` example since Gaussian has a lot more options the user might want to consider, but it's otherwise the same.

Like always, we have a `@job` decorator around a `make` function that defines the job. Here, however, note that there are optional `charge` and `mult` keywoard arguments in addition to the usual `self` and `atoms` arguments. Again, this is mostly for user convenience since it's something that is often specified.

Since there are a lot of arguments to specify, we define several defaults and then merge them with the user-specified arguments. This is done using the `merge_dicts` function, which is a helper function that is defined in `quacc.utils`. This function takes two dictionaries and merges them together, with the second dictionary taking precedence. It also has an option to remove any keys that have a value of `None` (which is what we do here).

Finally, we define the calculator, attach it to the `Atoms` object, run the calculation, and summarize the run. Here, we are not using ASE's optimizer routines so can use the regular `run_calc` and `summarize_run` functions in QuAcc. The optimization is selected here in the `defaults` dictionary with the `"opt"` key, as described the `Gaussian` calculator documentation.

This example nicely highlights how you'll often want to define some default values for a given job but also provide the user with full flexibility to change them as they see fit.
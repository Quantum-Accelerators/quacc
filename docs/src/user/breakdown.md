# Understanding the Code

## Pre-Requisites

The concepts in this section are based around [Jobflow](https://materialsproject.github.io/jobflow/), namely the concepts of a `job`, `flow`, `Maker`, and `Response`. It may be helpful to review the Jobflow tutorials if you are looking for more detailed information about any of these asepcts used in the examples below.

## Worked Examples

Before writing your own workflows, let's first break down a few examples already provided by QuAcc so that you understand the nuts and bolts.

### Example 1

For this example, I'll walk you through [the code](https://github.com/arosen93/quacc/blob/main/quacc/recipes/gaussian/core.py#L96-L164) for `quacc.recipes.tblite.core.StaticJob` step-by-step, which you should first skim over an try to understand what is going on.

Alright, now to explain a bit. Following the imports, the first notable lines are as follows:

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

Here, we have defined a job with the name `StaticJob`. Strictly speaking, this what is called a jobflow `Maker` class, which is a class that can be used to generate a `job`. You always start by inheriting from the `Maker` class in Jobflow. The `@dataclass` decorator allows us to define a class with attributes that are automatically initialized when the class is instantiated and should also always be included.

Below the name of the class, we have some docstrings to highlight the inputs that will define the settings of the job to be run. The first argument must always be `name` and is the name of the job. This is the only argument that needs to be included when writing a new `Maker`, but you can include as many additional arguments as you like. Here, we also include an optional argument `method` that will pass to the `tblite` calculator what level of theory to run. Finally, we have an optional argument `tblite_kwargs` that allows the user to pass a dictionary of any additional keyword arguments to the `tblite` calculator. It is important to give the user flexibility, so it is always a good idea to include a dictionary of arguments the user can pass on to the calculator if they so choose.

In case you are not familiar with the syntax here, `name: str = "tblite-Static"` means that the `name` argument takes an argument in the form of a string and by default, this is set to `"tblite-Static"`. The line `tblite_kwargs: Dict[str, Any] = field(default_factory=dict)` means that the `tblite_kwargs` optional argument takes a dictionary of any type of value, and by default, this is set to an empty dictionary. Of course, if you are coding this up only for internal use, you can hard-code the default values for the various arguments if you'd like, but this is good programming practice.

Okay, let's look at the latter half of the codeblock now:

```python
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

It's time to talk about the all-important `make` function. This is the function that will define the job itself using the aforementioned arguments. The `make` function must always take at least the following arguments: `self` and `atoms`. The former allows the job to inherit the arguments we defined above (i.e. `name`, `method`, `tblite_kwargs`). The latter is necessary so the job knows what `Atoms` object to act on. Again, if you aren't familiar with the syntax here, `atoms: Atoms` means that it is a required positional argument with type `Atoms` and has no default value (i.e. the user will supply it). Additionally, `-> Dict[str, Any]` is a typehint that tells the user what the output type will be. In this case, it is a dictionary of strings and any type of value. Typehints don't change the operation of the code; it is just helpful for the user to know what to expect.

The `@job` decorator is incredibly important. It tells Jobflow that the `make` function is a `job` and that it should be executed at runtime. All `job` objects need to have a `@job` wrapper preceeding it. The `@requires` decorator, on the other hand, is never absolutely needed. Here, it is included as a guide for the user; it will throw an error if `tblite` can't be imported, which typically implies it wasn't installed.

Below `make`, we have some docstrings again highlighting both the inputs and the outputs. Then we do four things:

- We instantiate the `TBLite` [ASE calculator](https://tblite.readthedocs.io/en/latest/users/ase.html) with the `method` and any `tblite_kwargs` arguments we defined above. The `**` before th `self.tblite_kwargs` is necessary to unpack the dictionary into individual keyword arguments, if there are any. We also attached this calculator to the `Atoms` object supplied by the user.

- We then run the calculation using the `run_calc` function supplied with QuAcc. Internally, this function is a wrapper around ASE's `.get_potential_energy()` function that will effectively launch a calculation given an `Atoms` object with an attached calculator. Crucially, by design, QuAcc does *not* modify the `Atoms` objects in place whenever possible, so we need to save the new `Atoms` object with a new variable name (here, `new_atoms` instead of `atoms`).

- We create a summary of the run inputs and outputs, which is a dictionary that will be output and stored. This is done using the `summarize_run` function supplied with QuAcc, which takes the output `Atoms`, input `Atoms`, and any other additional fields you'd like to specify. It is always recommended to store the `name` of the job as a field in the summary.

- We return the summary dictionary to the user.

That's it! Again, from the user perspective you would simplu create this job as `StaticJob().make(atoms)` once imported.

### Example 2

Let's continue with the `tblite` example!

This time, take a look at [the code](https://github.com/arosen93/quacc/blob/main/quacc/recipes/tblite/core.py#L69-L131) for `quacc.recipes.tblite.core.RelaxJob`, which is a geometry optimization.

What's different about this one? Well, the main thing is that some ASE calculators (including `tblite`) require you to use an ASE-based optimizer to relax the structure (see [here](https://wiki.fysik.dtu.dk/ase/ase/optimize.html) for more details). Just like how `run_calc` was a wrapper around ASE's `.get_potential_energy()` function to run a calculation, the `run_ase_opt` function is a wrapper around ASE's `Optimize` class. This makes it possible to run a geometry optimization in QuAcc using ASE. Otherwise, everything is largely the same except for some different arguments and a different summary function tailored for ASE optimization-based jobs.

### Example 3

Finally, let's look at a more typical relaxation job where you can run the relaxation using the code itself rather than ASE. For this example, take a look at [the code](https://github.com/arosen93/quacc/blob/main/quacc/recipes/gaussian/core.py#L96-L164) for `quacc.rescipes.gaussian.core.RelaxJob`, which is based around the [ASE Gaussian calculator](https://wiki.fysik.dtu.dk/ase/ase/calculators/gaussian.html). I explain the bits and pieces below.

Again, this should largely look similar. We start, like always, by instantiating a class that inherits a jobflow `Maker`. We then define the parameters that we want to be able to specify when we run the job. There are a few more parameters here than in the `tblite` example since Gaussian has a lot more options that we might want to suggest to the user that they should consider, but it's otherwise basically the same.

We have a `@job` decorator around a `make` function that defines the job. Here, however, note that there are optional `charge` and `mult` keywoard arguments in addition to the usual `self` and `atoms` arguments. Again, this is mostly for user convenience since it's something that is often specified in Gaussian calculations.

Since there are a lot of arguments to specify in a typical Gaussian relaxation calculation, we define several defaults and then merge them with any user-specified arguments. This is done using the `merge_dicts` function, which is a helper function that is defined in `quacc.utils`. This function takes two dictionaries and merges them together, with the second dictionary taking precedence. It also has an option to remove any keys that have a value of `None`, which is what we do here.

Finally, we define the calculator, attach it to the `Atoms` object, run the calculation, and summarize the run. Here, we are not using ASE's optimizer routines. Therefore, we can use the regular `run_calc` and `summarize_run` functions in QuAcc. The optimization is selected here in the `defaults` dictionary via the `"opt"` key, as described the `Gaussian` calculator documentation.

Importantly, this example highlights how you'll often want to define some default values for a given job but also provide the user with full flexibility to change them as they see fit. When designing new recipes, this is an important consideration that generally must be made to ensure that the recipe is as flexible as possible while still providing reasonable defaults.

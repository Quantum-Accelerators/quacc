
# Writing New Jobs and Flows

## Running Locally

The first step is to make a directory where you'll store your custom recipes. Let's make a folder that we'll name, for demonstration purposes, `myquacc`. Add a Python file to the `myquacc` directory that will hold your custom recipes, which we'll call `myrecipes.py`. We'll be editing `myrecipes.py` in the following examples.

For a custom job named `MyBeautifulJob`, this will mean that when you're done, you can import it as `from myquacc.myrecipes import MyBeautifulJob` provided the `myquacc` folder is either in your current working directory is added to your `PYTHONPATH` environment variable.

## Running with FireWorks or Contributing to QuAcc

If you are using FireWorks to manage your high-throughput workflows, all functions and classes to be installed in your environment in order for them to be successfully accessed at runtime. Therefore, instead of relying on a specific directory to store your custom recipes, you will need to add them to a copy of the QuAcc package itself. This is also necessary if you want to contribute new recipes to the main QuAcc repository.

In general, take the following steps:

1. Download the QuAcc repository from GitHub (e.g. `git clone https://github.com/arosen93/quacc.git`).
2. From the command line, change directories to the base directory of the QuAcc package you unpacked (e.g. `cd quacc`).
3. In your desired Python environment, install QuAcc in editable mode (i.e. `pip install -e .[all]`). This ensures that any changes you make to QuAcc will be instantly reflected in your code without having to uninstall/reinstall QuAcc.
4. Add your own recipes to the `quacc/recipes` directory, either as part of an existing folder or a new one.

If you added a new job (e.g. `MyBeautifulJob`) to a custom folder called `myrecipes`, then in your scripts you'd simply import it as `from quacc.recipes.myrecipes import MyBeautifulJob`, similar to what you might do with the pre-packaged recipes.

## Adding a Job or Flow for a New Calculator

While QuAcc ships with support for a number of common quantum chemistry packages, there may be others you wish to use. Natively, QuAcc can easily support can calculator that is part of the [ASE ecosystem](https://wiki.fysik.dtu.dk/ase/ase/calculators/calculators.html), of which there are many! If you want to use a calculator that is not part of the ASE ecosystem, you can still use QuAcc, but you'll need to write an ASE calculator for it first (or at least a Python API of some kind).

## Generating an Output Schema

At the end of any recipe you make, you need to return a summary dictionary of the main inputs and outputs. As alluded to in [breakdown.md](breakdown.md), this can currently be done with one of several possible functions:

1. `quacc.schemas.calc.summarize_run`: This function is used to summarize any typical run that is done through ASE. It summarizes the inputs and outputs from the ASE calculator object. This is the most common function to use.

2. `quacc.schemas.cclib.summarize_run`: This function is used to summarize the inputs and outputs of a run for a code that is supported by [cclib](https://cclib.github.io/data.html). This is preferred over `quacc.schemas.calc.summarize_run` if the code is supported by cclib because it will generate a more detailed summary than ASE can provide.

3. `quacc.schemas.calc.summarize_opt_run`: This function is used to summarize any relaxation run that is specifically done via an ASE `Optimizer`. This is the only way ASE `Optimizer` runs can be summarized.

4. `quacc.schemas.vasp.summarize_run`: This function is used to summarize a VASP run. It has been custom-made for VASP.

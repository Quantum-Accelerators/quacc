
# Writing New Jobs and Flows

Oftentimes, you will want to write your own jobs and flows for QuAcc, whether they are for your own internal use or to contribute back to the QuAcc package. This page will walk you through the process of writing your own jobs and flows.

## Method 1: Installing QuAcc Recipes

Perhaps the preferred way to write your recipes in QuAcc is to install the code in editable mode and add any new recipes to `quacc.recipes`. This method is described below. This approach is particularly necessary if you plan to contribute your recipes to the main QuAcc repository. It is also necessary if you are using FireWorks to manage your high-throughput workflows, as all functions and classes must be installed in your environment in order for them to be successfully accessed at runtime.

In general, take the following steps:

1. Download the QuAcc repository from GitHub (e.g. `git clone https://github.com/arosen93/quacc.git`).
2. From the command line, change directories to the base directory of the QuAcc package you unpacked (e.g. `cd quacc`).
3. In your desired Python environment, install QuAcc in editable mode (i.e. `pip install -e .[all]`). This ensures that any changes you make to QuAcc will be instantly reflected in your code without having to uninstall/reinstall QuAcc.
4. Add your own recipes to the `quacc/recipes` directory, either as part of an existing folder or a new one.

If you added a new job (e.g. `MyBeautifulJob`) to a custom folder called `myrecipes`, then in your scripts you'd simply import it as `from quacc.recipes.myrecipes import MyBeautifulJob`, similar to what you might do with the pre-packaged recipes.

## Method 2: Storing QuAcc Recipes Locally

Alternatively, if you don't want to install QuAcc in editable mode, you can still write your own recipes and use them as you like, provided you are not using FireWorks.

Instead of storing the recipes in the QuAcc code itself, the following steps should be carried out:

1. Make a custom directory where you'll store your custom recipes. For instance, you could make a folder called `myquacc`.
2. Add a Python file to the `myquacc` directory that will hold your custom recipes, which you could call `myrecipes.py`.
3. Add your custom recipes to the `myrecipes.py` file, drawing paralells with how they are constructed in `quacc.recipes` for inspiration. Let's pretend you made a new job called `MyBeautifulJob`.
4. When you're done, you can call this newly made recipe from any Python script by importing it as `from myquacc.myrecipes import MyBeautifulJob`, provided the `myquacc` folder is either in your current working directory or is added to your `PYTHONPATH` environment variable.

## Adding a Job or Flow for a New Calculator

While QuAcc ships with support for a number of common quantum chemistry and computational materials science packages, there may be others you wish to use. Natively, QuAcc can easily support can calculator that is part of the [ASE ecosystem](https://wiki.fysik.dtu.dk/ase/ase/calculators/calculators.html), of which there are many to pick from. If you want to use a calculator that is not part of the ASE ecosystem, you can still use QuAcc, but you'll need to write an ASE calculator for it first (or at least a Python API of some kind).

## Generating an Output Schema

At the end of any recipe you make, you need to return a summary dictionary of the main inputs and outputs. As alluded to in ["Understanding the Code"](breakdown.md), this can currently be done with one of several possible functions:

1. `quacc.schemas.calc.summarize_run`: This function is used to summarize a typical run that is done through ASE. It summarizes the inputs and outputs by pulling them from the ASE calculator and `Atoms` object. This is the most common function to use.

2. `quacc.schemas.cclib.summarize_run`: This function is used to summarize the inputs and outputs of a run for codes that are supported by [cclib](https://cclib.github.io/data.html). When possible, this is preferred over `quacc.schemas.calc.summarize_run` because cclib will generate a more detailed summary than ASE can provide. It is also more standardized across different codes.

3. `quacc.schemas.calc.summarize_opt_run`: This function is used to summarize a relaxation run that is specifically carried out via an ASE `Optimizer`. This is the only way ASE `Optimizer` runs can be summarized within QuAcc.

4. `quacc.schemas.vasp.summarize_run`: This function is used to summarize a VASP run. It has been custom-made for VASP based on the schema used in `Atomate2`.

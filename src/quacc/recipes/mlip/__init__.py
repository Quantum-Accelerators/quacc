"""
Recipes for machine-learned interatomic potentials.

!!! Warning

    The MLIP calculators are cached and shared by all jobs in the same Python
    process, so MLIP jobs must not run concurrently in multiple threads of one
    process. See [quacc.recipes.mlip._base.pick_calculator][] for details.
"""

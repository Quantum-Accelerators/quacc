# Connecting Quacc with Atomate2

[Atomate2](https://github.com/materialsproject/atomate2) is a fantastic computational materials science workflow program that shares many similarities with quacc. If you wish to combine workflows from quacc with those from Atomate2, that is possible through the use of Jobflow.

!!! Tip

    Refer to the [Using a Workflow Engine](wflow_engines.md) guide for more information on how to use Jobflow with quacc.

For instance, consider a toy example where you want to relax a bulk Cu structure using the `tblite` calculator in quacc, which you then use to run a VASP relaxation and bandstructure workflow in Atomate2. That would look like the following.

```python
from ase.build import bulk
from atomate2.vasp.flows.core import RelaxBandStructureMaker
from jobflow import Flow
from quacc.recipes.tblite.core import relax_job

atoms = bulk("Cu")

job1 = relax_job(atoms)
bandstructure_flow = RelaxBandStructureMaker().make_flow(job1.output["structure"]) # (1)!

flow = Flow([job1]) + bandstructure_flow # (2)!
```

1. All Atomate2 workflows take a Pymatgen `Structure` or `Molecule` object as input. This is one of the properties in the returned output of a quacc recipe, which is why we can do `job1.output["structure"]`.

2. The `+` operator can be used to combine two flows into one. We converted the first job into its own `Flow` definition to enable this.

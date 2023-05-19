import os
from shutil import rmtree

import pytest
from ase.build import bulk, molecule

from quacc.recipes.emt.core import relax_job, static_job
from quacc.recipes.emt.slabs import BulkToSlabsFlow

try:
    import covalent as ct
except ImportError:
    ct = None


def teardown_module():
    for f in os.listdir(os.getcwd()):
        if (
            f.endswith(".log")
            or f.endswith(".pckl")
            or f.endswith(".traj")
            or f.endswith(".out")
            or ".gz" in f
        ):
            os.remove(f)
        if "quacc-tmp" in f or f == "tmp_dir":
            rmtree(f)


@pytest.mark.skipif(
    ct is None,
    reason="Covalent must be installed",
)
def test_covalent_config():
    ct_config = ct.get_config()
    for executor in ["local", "dask"]:
        assert (
            ct_config["executors"][executor].get("create_unique_workdir", False) is True
        )
    if "slurm" in ct_config["executors"]:
        assert ct_config["executors"]["slurm"].get("use_srun", True) is False
        assert (
            ct_config["executors"]["slurm"].get("create_unique_workdir", False) is True
        )


@pytest.mark.skipif(
    ct is None or os.environ.get("GITHUB_ACTIONS", False) is False,
    reason="This test is only meant to be run on GitHub Actions",
)
def test_tutorials():
    # Test of the various tutorials

    @ct.lattice
    def workflow(atoms):
        result1 = relax_job(atoms)
        result2 = static_job(result1["atoms"])
        return result2

    atoms = bulk("Cu")
    dispatch_id = ct.dispatch(workflow)(atoms)
    result = ct.get_result(dispatch_id, wait=True)
    assert result.status == "COMPLETED"

    # ------------------------------------------------------------
    @ct.lattice
    def workflow1(atoms1, atoms2):
        result1 = relax_job(atoms1)
        result2 = relax_job(atoms2)

        return [result1, result2]

    atoms1 = bulk("Cu")
    atoms2 = molecule("N2")
    dispatch_id = ct.dispatch(workflow1)(atoms1, atoms2)
    result = ct.get_result(dispatch_id, wait=True)
    assert result.status == "COMPLETED"

    # ------------------------------------------------------------
    @ct.lattice
    def workflow2(atoms):
        relaxed_bulk = relax_job(atoms)
        relaxed_slabs = BulkToSlabsFlow(slab_static_electron=None).run(
            relaxed_bulk["atoms"]
        )
        return relaxed_slabs

    atoms = bulk("Cu")
    dispatch_id = ct.dispatch(workflow2)(atoms)
    result = ct.get_result(dispatch_id, wait=True)
    assert result.status == "COMPLETED"

    # ------------------------------------------------------------
    @ct.lattice(executor="local")
    def workflow3(atoms):
        result1 = relax_job(atoms)
        result2 = static_job(result1["atoms"])
        return result2

    atoms = bulk("Cu")
    dispatch_id = ct.dispatch(workflow3)(atoms)
    result = ct.get_result(dispatch_id, wait=True)
    assert result.status == "COMPLETED"

    # ------------------------------------------------------------
    @ct.electron(executor="dask")
    def relax_electron(atoms):
        return relax_job(atoms)

    @ct.electron(executor="local")
    def static_electron(atoms):
        return static_job(atoms)

    @ct.lattice
    def workflow4(atoms):
        output1 = relax_electron(atoms)
        output2 = static_electron(output1["atoms"])
        return output2

    atoms = bulk("Cu")
    dispatch_id = ct.dispatch(workflow4)(atoms)
    result = ct.get_result(dispatch_id, wait=True)
    assert result.status == "COMPLETED"

    # ------------------------------------------------------------
    @ct.lattice
    def workflow5(atoms):
        job1 = relax_job
        job1.electron_object.metadata["executor"] = "dask"

        job2 = static_job
        job2.electron_object.metadata["executor"] = "local"

        output1 = job1(atoms)
        output2 = job2(output1["atoms"])
        return output2

    atoms = bulk("Cu")
    dispatch_id = ct.dispatch(workflow5)(atoms)
    result = ct.get_result(dispatch_id, wait=True)
    assert result.status == "COMPLETED"

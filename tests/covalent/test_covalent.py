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
    assert ct_config["executors"]["local"].get("create_unique_workdir", False) is True
    if "dask" in ct_config["executors"]:
        assert (
            ct_config["executors"]["dask"].get("create_unique_workdir", False) is True
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
        job1 = ct.electron(relax_job)
        result1 = job1(atoms)
        job2 = ct.electron(static_job)
        result2 = job2(result1["atoms"])
        return result2

    atoms = bulk("Cu")
    dispatch_id = ct.dispatch(workflow)(atoms)
    result = ct.get_result(dispatch_id, wait=True)
    assert result.status == "COMPLETED"

    # ------------------------------------------------------------
    @ct.lattice
    def workflow(atoms1, atoms2):
        result1 = ct.electron(relax_job)(atoms1)
        result2 = ct.electron(relax_job)(atoms2)

        return [result1, result2]

    atoms1 = bulk("Cu")
    atoms2 = molecule("N2")
    dispatch_id = ct.dispatch(workflow)(atoms1, atoms2)
    result = ct.get_result(dispatch_id, wait=True)
    assert result.status == "COMPLETED"

    # ------------------------------------------------------------
    # @ct.lattice
    # def workflow(atoms):
    #     relaxed_bulk = ct.electron(relax_job)(atoms)
    #     relaxed_slabs = BulkToSlabsFlow(static_electron=None).run(relaxed_bulk["atoms"])
    #     return relaxed_slabs

    # atoms = bulk("Cu")
    # dispatch_id = ct.dispatch(workflow)(atoms)
    # result = ct.get_result(dispatch_id, wait=True)
    # assert result.status == "COMPLETED"

    # ------------------------------------------------------------
    @ct.lattice(executor="local")
    def workflow(atoms):
        result1 = ct.electron(relax_job)(atoms)
        result2 = ct.electron(static_job)(result1["atoms"])
        return result2

    atoms = bulk("Cu")
    dispatch_id = ct.dispatch(workflow)(atoms)
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
    def workflow(atoms):
        output1 = relax_electron(atoms)
        output2 = static_electron(output1["atoms"])
        return output2

    atoms = bulk("Cu")
    dispatch_id = ct.dispatch(workflow)(atoms)
    result = ct.get_result(dispatch_id, wait=True)
    assert result.status == "COMPLETED"

    # ------------------------------------------------------------
    @ct.lattice
    def workflow(atoms):
        job1 = ct.electron(relax_job)
        job1.electron_object.metadata["executor"] = "dask"

        job2 = ct.electron(static_job)
        job2.electron_object.metadata["executor"] = "local"

        output1 = job1(atoms)
        output2 = job2(output1["atoms"])
        return output2

    atoms = bulk("Cu")
    dispatch_id = ct.dispatch(workflow)(atoms)
    result = ct.get_result(dispatch_id, wait=True)
    assert result.status == "COMPLETED"

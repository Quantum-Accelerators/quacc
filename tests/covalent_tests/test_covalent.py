import pytest
from ase.build import bulk

from quacc import SETTINGS, flow, job, subflow
from quacc.recipes.emt.core import relax_job, static_job

try:
    import covalent as ct

    ct = ct if SETTINGS.WORKFLOW_ENGINE == "covalent" else None
except ImportError:
    ct = None


@pytest.mark.skipif(ct is None, reason="This test requires Covalent")
def test_covalent_decorators(tmpdir):
    tmpdir.chdir()

    @job
    def add(a, b):
        return a + b

    @job
    def mult(a, b):
        return a * b

    @job
    def make_more(val):
        return [val] * 3

    @subflow
    def add_distributed(vals, c):
        return [add(val, c) for val in vals]

    @flow
    def workflow(a, b, c):
        return mult(add(a, b), c)

    @flow
    def dynamic_workflow(a, b, c):
        result1 = add(a, b)
        result2 = make_more(result1)
        return add_distributed(result2, c)

    assert add(1, 2) == 3
    assert mult(1, 2) == 2
    assert ct.get_result(ct.dispatch(workflow)(1, 2, 3), wait=True).result == 9
    assert ct.get_result(ct.dispatch(dynamic_workflow)(1, 2, 3), wait=True).result == [
        6,
        6,
        6,
    ]
    assert ct.get_result(
        ct.dispatch(flow(add_distributed))([1, 1, 1], 2), wait=True
    ).result == [
        3,
        3,
        3,
    ]


@pytest.mark.skipif(ct is None, reason="This test requires Covalent")
def test_covalent_decorators_args(tmpdir):
    tmpdir.chdir()

    @job(executor="local")
    def add(a, b):
        return a + b

    @job(executor="local")
    def mult(a, b):
        return a * b

    @job(executor="local")
    def make_more(val):
        return [val] * 3

    @subflow(executor="local")
    def add_distributed(vals, c):
        return [add(val, c) for val in vals]

    @flow(executor="local")
    def workflow(a, b, c):
        return mult(add(a, b), c)

    @flow(executor="local")
    def dynamic_workflow(a, b, c):
        result1 = add(a, b)
        result2 = make_more(result1)
        return add_distributed(result2, c)

    assert add(1, 2) == 3
    assert mult(1, 2) == 2
    assert ct.get_result(ct.dispatch(workflow)(1, 2, 3), wait=True).result == 9
    assert ct.get_result(ct.dispatch(dynamic_workflow)(1, 2, 3), wait=True).result == [
        6,
        6,
        6,
    ]


@pytest.mark.skipif(ct is None, reason="This test requires Covalent")
def test_quickstart(tmpdir):
    tmpdir.chdir()

    import covalent as ct
    from ase.build import bulk

    from quacc.recipes.emt.slabs import bulk_to_slabs_flow

    # Define the Atoms object
    atoms = bulk("Cu")

    # Dispatch the workflow
    dispatch_id = ct.dispatch(bulk_to_slabs_flow)(atoms)

    # Fetch the results
    result = ct.get_result(dispatch_id, wait=True)
    assert result.status == "COMPLETED"


@pytest.mark.skipif(ct is None, reason="This test requires Covalent")
def test_tutorial1a(tmpdir):
    tmpdir.chdir()

    import covalent as ct
    from ase.build import bulk

    from quacc import flow
    from quacc.recipes.emt.core import relax_job

    # Make an Atoms object of a bulk Cu structure
    atoms = bulk("Cu")

    # Define the workflow
    workflow = flow(relax_job)  # (1)!

    # Dispatch the workflow to the Covalent server
    # with the bulk Cu Atoms object as the input
    dispatch_id = ct.dispatch(workflow)(atoms)  # (3)!

    # Fetch the result from the server
    result = ct.get_result(dispatch_id, wait=True)  # (4)!
    assert result.status == "COMPLETED"


@pytest.mark.skipif(ct is None, reason="This test requires Covalent")
def test_tutorial1b(tmpdir):
    tmpdir.chdir()

    import covalent as ct
    from ase.build import bulk

    from quacc.recipes.emt.slabs import bulk_to_slabs_flow

    atoms = bulk("Cu")
    dispatch_id = ct.dispatch(bulk_to_slabs_flow)(atoms)  # (1)!
    result = ct.get_result(dispatch_id, wait=True)
    assert result.status == "COMPLETED"


@pytest.mark.skipif(ct is None, reason="This test requires Covalent")
def test_tutorial2a(tmpdir):
    tmpdir.chdir()
    import covalent as ct
    from ase.build import bulk

    from quacc import flow
    from quacc.recipes.emt.core import relax_job, static_job

    # Define the workflow
    @flow  # (1)!
    def workflow(atoms):
        # Define Job 1
        result1 = relax_job(atoms)  # (2)!

        # Define Job 2, which takes the output of Job 1 as input
        return static_job(result1)

    # Make an Atoms object of a bulk Cu structure
    atoms = bulk("Cu")

    # Dispatch the workflow to the Covalent server
    # with the bulk Cu Atoms object as the input
    dispatch_id = ct.dispatch(workflow)(atoms)  # (3)!

    # Fetch the result from the server
    result = ct.get_result(dispatch_id, wait=True)  # (4)!
    assert result.status == "COMPLETED"


@pytest.mark.skipif(ct is None, reason="This test requires Covalent")
def test_tutorial2b(tmpdir):
    tmpdir.chdir()

    from ase.build import bulk, molecule

    from quacc import flow
    from quacc.recipes.emt.core import relax_job

    # Define workflow
    @flow
    def workflow(atoms1, atoms2):
        # Define two independent relaxation jobs
        result1 = relax_job(atoms1)
        result2 = relax_job(atoms2)

        return {"result1": result1, "result2": result2}

    # Define two Atoms objects
    atoms1 = bulk("Cu")
    atoms2 = molecule("N2")

    # Dispatch the workflow to the Covalent server
    dispatch_id = ct.dispatch(workflow)(atoms1, atoms2)

    # Fetch the results from the server
    result = ct.get_result(dispatch_id, wait=True)
    assert result.status == "COMPLETED"


@pytest.mark.skipif(ct is None, reason="This test requires Covalent")
def test_tutorial2c(tmpdir):
    tmpdir.chdir()

    import covalent as ct
    from ase.build import bulk

    from quacc import flow
    from quacc.recipes.emt.core import relax_job
    from quacc.recipes.emt.slabs import bulk_to_slabs_flow

    @flow
    def workflow(atoms):
        relaxed_bulk = relax_job(atoms)
        return bulk_to_slabs_flow(relaxed_bulk, run_static=False)  # (1)!

    atoms = bulk("Cu")
    dispatch_id = ct.dispatch(workflow)(atoms)
    result = ct.get_result(dispatch_id, wait=True)
    assert result.status == "COMPLETED"


@pytest.mark.skipif(ct is None, reason="This test requires Covalent")
def test_tutorial_excecutor1(tmpdir):
    tmpdir.chdir()

    @flow(executor="local")
    def workflow4(atoms):
        result1 = relax_job(atoms)
        return static_job(result1)

    atoms = bulk("Cu")
    dispatch_id = ct.dispatch(workflow4)(atoms)
    result = ct.get_result(dispatch_id, wait=True)
    assert result.status == "COMPLETED"


@pytest.mark.skipif(ct is None, reason="This test requires Covalent")
def test_tutorial_excecutor2(tmpdir):
    tmpdir.chdir()

    relax_job.electron_object.executor = "dask"
    static_job.electron_object.executor = "local"

    @flow
    def workflow5(atoms):
        output1 = relax_job(atoms)
        return static_job(output1)

    atoms = bulk("Cu")
    dispatch_id = ct.dispatch(workflow5)(atoms)
    result = ct.get_result(dispatch_id, wait=True)
    assert result.status == "COMPLETED"


@pytest.mark.skipif(ct is None, reason="This test requires Covalent")
def test_comparison1(tmpdir):
    tmpdir.chdir()

    import covalent as ct

    from quacc import flow, job

    @job  # (1)!
    def add(a, b):
        return a + b

    @job
    def mult(a, b):
        return a * b

    @flow  # (2)!
    def workflow(a, b, c):
        return mult(add(a, b), c)

    # Locally
    result = workflow(1, 2, 3)  # 9  (3)!

    # Dispatched
    dispatch_id = ct.dispatch(workflow)(1, 2, 3)  # (4)!
    result = ct.get_result(dispatch_id, wait=True)  # 9  (5)!
    assert result.status == "COMPLETED"


@pytest.mark.skipif(ct is None, reason="This test requires Covalent")
def test_comparison2(tmpdir):
    tmpdir.chdir()

    import covalent as ct

    from quacc import flow, job, subflow

    @job
    def add(a, b):
        return a + b

    @job
    def make_more(val):
        return [val] * 3

    @subflow  # (1)!
    def add_distributed(vals, c):
        return [add(val, c) for val in vals]

    @flow
    def workflow(a, b, c):
        result1 = add(a, b)
        result2 = make_more(result1)
        return add_distributed(result2, c)

    # Locally
    result = workflow(1, 2, 3)  # e.g. [6, 6, 6]

    # Dispatched
    dispatch_id = ct.dispatch(workflow)(1, 2, 3)
    result = ct.get_result(dispatch_id, wait=True)  # e.g. [6, 6, 6]

    assert result.status == "COMPLETED"


@pytest.mark.skipif(ct is None, reason="This test requires Covalent")
def test_comparison3(tmpdir):
    tmpdir.chdir()
    import covalent as ct

    from quacc import flow, job

    @job  #  (1)!
    def add(a, b):
        return a + b

    @job
    def mult(a, b):
        return a * b

    @flow  #  (2)!
    def workflow(a, b, c):
        return mult(add(a, b), c)

    dispatch_id = ct.dispatch(workflow)(1, 2, 3)  # (3)!
    result = ct.get_result(dispatch_id, wait=True)
    assert result.status == "COMPLETED"


@pytest.mark.skipif(ct is None, reason="This test requires Covalent")
def test_comparison4(tmpdir):
    tmpdir.chdir()

    import covalent as ct

    from quacc import flow, job, subflow

    @job
    def add(a, b):
        return a + b

    @job
    def make_more(val):
        return [val] * 3

    @subflow  #  (1)!
    def add_distributed(vals, c):
        return [add(val, c) for val in vals]

    @flow
    def workflow(a, b, c):
        result1 = add(a, b)
        result2 = make_more(result1)
        return add_distributed(result2, c)

    # Dispatched
    dispatch_id = ct.dispatch(workflow)(1, 2, 3)
    result = ct.get_result(dispatch_id, wait=True)  # e.g. [6, 6, 6]

    assert result.status == "COMPLETED"

import contextlib

import pytest

from quacc import SETTINGS, flow, job, subflow

try:
    import parsl

    parsl = parsl if SETTINGS.WORKFLOW_ENGINE == "parsl" else None

except ImportError:
    parsl = None
try:
    import psi4
except ImportError:
    psi4 = None
try:
    from tblite.ase import TBLite
except ImportError:
    TBLite = None
try:
    from quacc.recipes.emt.defects import bulk_to_defects_flow
except ImportError:
    bulk_to_defects_flow = None
DEFAULT_SETTINGS = SETTINGS.copy()


def setup_module():
    if parsl:
        with contextlib.suppress(Exception):
            parsl.load()


@pytest.mark.skipif(parsl is None, reason="Parsl not installed")
def test_parsl_decorators(tmpdir):
    tmpdir.chdir()
    SETTINGS.WORKFLOW_ENGINE = "parsl"

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

    assert add(1, 2).result() == 3
    assert mult(1, 2).result() == 2
    assert workflow(1, 2, 3).result() == 9
    assert dynamic_workflow(1, 2, 3).result() == [6, 6, 6]


@pytest.mark.skipif(parsl is None, reason="Parsl not installed")
def test_parsl_decorators_args(tmpdir):
    tmpdir.chdir()
    SETTINGS.WORKFLOW_ENGINE = "parsl"

    @job()
    def add(a, b):
        return a + b

    @job()
    def mult(a, b):
        return a * b

    @job()
    def make_more(val):
        return [val] * 3

    @subflow()
    def add_distributed(vals, c):
        return [add(val, c) for val in vals]

    @flow()
    def workflow(a, b, c):
        return mult(add(a, b), c)

    @flow()
    def dynamic_workflow(a, b, c):
        result1 = add(a, b)
        result2 = make_more(result1)
        return add_distributed(result2, c)

    assert add(1, 2).result() == 3
    assert mult(1, 2).result() == 2
    assert workflow(1, 2, 3).result() == 9
    assert dynamic_workflow(1, 2, 3).result() == [6, 6, 6]


@pytest.mark.skipif(
    parsl is None,
    reason="Parsl is not installed or specified in config",
)
def test_tutorial1a(tmpdir):
    tmpdir.chdir()

    from ase.build import bulk

    from quacc.recipes.emt.core import relax_job

    # Make an Atoms object of a bulk Cu structure
    atoms = bulk("Cu")

    # Call the PythonApp
    future = relax_job(atoms)  # (1)!

    # Print result
    assert "atoms" in future.result()  # (2)!


@pytest.mark.skipif(
    parsl is None,
    reason="Parsl is not installed or specified in config",
)
def test_tutorial1b(tmpdir):
    tmpdir.chdir()

    from ase.build import bulk

    from quacc.recipes.emt.core import relax_job

    # Make an Atoms object of a bulk Cu structure
    atoms = bulk("Cu")

    # Call the PythonApp
    future = relax_job(atoms)  # (1)!

    # Print result
    assert "atoms" in future.result()  # (2)!

    from ase.build import bulk

    from quacc.recipes.emt.slabs import bulk_to_slabs_flow

    # Define the Atoms object
    atoms = bulk("Cu")

    # Define the workflow
    future = bulk_to_slabs_flow(atoms)  # (1)!

    # Print the results
    assert "atoms" in future.result()[0]  # (2)!


@pytest.mark.skipif(
    parsl is None,
    reason="Parsl is not installed or specified in config",
)
def test_tutorial2a(tmpdir):
    tmpdir.chdir()

    from ase.build import bulk

    from quacc.recipes.emt.core import relax_job, static_job

    # Define the workflow
    def workflow(atoms):
        # Define Job 1
        future1 = relax_job(atoms)  # (1)!

        # Define Job 2, which takes the output of Job 1 as input
        return static_job(future1)

    # Make an Atoms object of a bulk Cu structure
    atoms = bulk("Cu")

    # Dispatch the workflow
    future = workflow(atoms)

    # Fetch the result
    result = future.result()  # (2)!
    assert "atoms" in result


@pytest.mark.skipif(
    parsl is None,
    reason="Parsl is not installed or specified in config",
)
def test_tutorial2b(tmpdir):
    tmpdir.chdir()

    from ase.build import bulk, molecule

    from quacc.recipes.emt.core import relax_job

    # Define workflow
    def workflow(atoms1, atoms2):
        # Define two independent relaxation jobs
        result1 = relax_job(atoms1)
        result2 = relax_job(atoms2)

        return {"result1": result1, "result2": result2}

    # Define two Atoms objects
    atoms1 = bulk("Cu")
    atoms2 = molecule("N2")

    # Define two independent relaxation jobs
    futures = workflow(atoms1, atoms2)

    # Fetch the results
    result1 = futures["result1"].result()
    result2 = futures["result2"].result()

    # Print the results
    assert "atoms" in result1
    assert "atoms" in result2


@pytest.mark.skipif(
    parsl is None,
    reason="Parsl is not installed or specified in config",
)
def test_tutorial2c(tmpdir):
    tmpdir.chdir()

    from ase.build import bulk

    from quacc.recipes.emt.core import relax_job
    from quacc.recipes.emt.slabs import bulk_to_slabs_flow

    # Define the workflow
    def workflow(atoms):
        relaxed_bulk = relax_job(atoms)
        return bulk_to_slabs_flow(relaxed_bulk, run_static=False)  # (1)!

    # Define the Atoms object
    atoms = bulk("Cu")

    # Dispatch the workflow
    future = workflow(atoms)

    # Fetch the results
    result = future.result()

    # Print the results
    assert len(result) == 4


@pytest.mark.skipif(
    parsl is None,
    reason="Parsl is not installed or specified in config",
)
def test_comparison1(tmpdir):
    tmpdir.chdir()

    from quacc import job

    @job  #  (1)!
    def add(a, b):
        return a + b

    @job
    def mult(a, b):
        return a * b

    def workflow(a, b, c):  #  (2)!
        return mult(add(a, b), c)

    result = workflow(1, 2, 3).result()  # 9
    assert result == 9


@pytest.mark.skipif(
    parsl is None,
    reason="Parsl is not installed or specified in config",
)
def test_comparison2(tmpdir):
    tmpdir.chdir()

    from quacc import job, subflow

    @job
    def add(a, b):
        return a + b

    @job
    def make_more(val):
        return [val] * 3

    @subflow  # (1)!
    def add_distributed(vals, c):
        return [add(val, c) for val in vals]

    future1 = add(1, 2)
    future2 = make_more(future1)
    future3 = add_distributed(future2, 3)

    assert future3.result() == [6, 6, 6]


@pytest.mark.skipif(
    parsl is None,
    reason="Parsl is not installed or specified in config",
)
def test_comparison3(tmpdir):
    tmpdir.chdir()
    from quacc import job

    @job  #  (1)!
    def add(a, b):
        return a + b

    @job
    def mult(a, b):
        return a * b

    future1 = add(1, 2)
    future2 = mult(future1, 3)

    assert future2.result() == 9


@pytest.mark.skipif(
    parsl is None,
    reason="Parsl is not installed or specified in config",
)
def test_comparison4(tmpdir):
    tmpdir.chdir()
    from quacc import job, subflow

    @job
    def add(a, b):
        return a + b

    @job
    def make_more(val):
        return [val] * 3

    @subflow  #  (1)!
    def add_distributed(vals, c):
        return [add(val, c) for val in vals]

    future1 = add(1, 2)
    future2 = make_more(future1)
    future3 = add_distributed(future2, 3)

    assert future3.result() == [6, 6, 6]


@pytest.mark.skipif(
    parsl is None,
    reason="Parsl is not installed or specified in config",
)
def test_docs_recipes_emt(tmpdir):
    tmpdir.chdir()

    from ase.build import bulk

    from quacc.recipes.emt.core import static_job

    atoms = bulk("Cu")
    future = static_job(atoms)
    result = future.result()
    assert future.done()

    # -----------------

    from ase.build import bulk

    from quacc.recipes.emt.slabs import bulk_to_slabs_flow

    atoms = bulk("Ni")
    future = bulk_to_slabs_flow(atoms)
    result = future.result()
    assert future.done()


@pytest.mark.skipif(
    parsl is None or bulk_to_defects_flow is None,
    reason="Parsl is not installed or specified in config",
)
def test_docs_recipes_emt_defects(tmpdir):
    tmpdir.chdir()
    from ase.build import bulk

    from quacc.recipes.emt.defects import bulk_to_defects_flow

    atoms = bulk("Cu")
    future = bulk_to_defects_flow(atoms)
    result = future.result()
    assert future.done()


@pytest.mark.skipif(
    parsl is None,
    reason="Parsl is not installed or specified in config",
)
def test_docs_recipes_lj(tmpdir):
    tmpdir.chdir()
    from ase.build import molecule

    from quacc.recipes.lj.core import relax_job

    atoms = molecule("N2")
    future = relax_job(atoms)
    result = future.result()
    assert future.done()

    # ---------------------
    from ase.build import molecule

    from quacc.recipes.lj.core import static_job

    atoms = molecule("N2")
    future = static_job(atoms)
    result = future.result()
    assert future.done()

    # --------------------
    from ase.build import molecule

    from quacc.recipes.lj.core import freq_job

    atoms = molecule("N2")
    future = freq_job(atoms)
    result = future.result()
    assert future.done()


@pytest.mark.skipif(
    parsl is None or psi4 is None,
    reason="Parsl is not installed or specified in config",
)
def test_docs_recipes_psi4(tmpdir):
    tmpdir.chdir()
    from ase.build import molecule

    from quacc.recipes.psi4.core import static_job

    atoms = molecule("O2")
    future = static_job(atoms, 0, 3, method="wb97m-v", basis="def2-svp")
    result = future.result()
    assert result.done()


@pytest.mark.skipif(
    parsl is None or TBLite is None,
    reason="Parsl is not installed or specified in config",
)
def test_docs_recipes_tblite(tmpdir):
    tmpdir.chdir()

    from ase.build import bulk

    from quacc.recipes.tblite.core import relax_job

    atoms = bulk("C")
    future = relax_job(atoms, relax_cell=True)
    result = future.result()
    assert future.done()

    # ------------------------

    from ase.build import bulk

    from quacc.recipes.tblite.core import static_job

    atoms = bulk("C")
    future = static_job(atoms, method="GFN1-xTB")
    result = future.result()
    assert future.done()

    # ------------------------
    from ase.build import molecule

    from quacc.recipes.tblite.core import freq_job

    atoms = molecule("N2")
    future = freq_job(atoms)
    result = future.result()
    assert future.done()

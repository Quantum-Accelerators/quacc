import pytest

from quacc import SETTINGS

try:
    import parsl
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
DEFAULT_SETTINGS = SETTINGS.copy()


def setup_module():
    if parsl:
        try:
            parsl.load()
        except RuntimeError:
            pass

    SETTINGS.WORKFLOW_ENGINE = "parsl"


def teardown_module():
    SETTINGS.WORKFLOW_ENGINE = DEFAULT_SETTINGS.WORKFLOW_ENGINE


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

    # Make an Atoms object of a bulk Cu structure
    atoms = bulk("Cu")

    # Call App 1
    future1 = relax_job(atoms)  # (1)!

    # Call App 2, which takes the output of App 1 as input
    future2 = static_job(future1)

    # Print result
    assert "atoms" in future2.result()


@pytest.mark.skipif(
    parsl is None,
    reason="Parsl is not installed or specified in config",
)
def test_tutorial2b(tmpdir):
    tmpdir.chdir()

    from ase.build import bulk, molecule

    from quacc.recipes.emt.core import relax_job

    # Define two Atoms objects
    atoms1 = bulk("Cu")
    atoms2 = molecule("N2")

    # Define two independent relaxation jobs
    future1 = relax_job(atoms1)
    future2 = relax_job(atoms2)

    # Print the results
    assert "atoms" in future1.result()
    assert "atoms" in future2.result()


@pytest.mark.skipif(
    parsl is None,
    reason="Parsl is not installed or specified in config",
)
def test_tutorial2c(tmpdir):
    tmpdir.chdir()
    from ase.build import bulk

    from quacc.recipes.emt.core import relax_job
    from quacc.recipes.emt.slabs import bulk_to_slabs_flow

    # Define the Atoms object
    atoms = bulk("Cu")

    # Define the workflow
    future1 = relax_job(atoms)
    future2 = bulk_to_slabs_flow(future1, slab_static=None)  # (1)!

    # Print the results
    assert len(future2.result()) == 4
    assert future2.done()


@pytest.mark.skipif(
    parsl is None,
    reason="Parsl is not installed or specified in config",
)
def test_comparison1(tmpdir):
    tmpdir.chdir()

    from quacc import job

    @job  # (1)!
    def add(a, b):
        return a + b

    @job
    def mult(a, b):
        return a * b

    future1 = add(1, 2)
    future2 = mult(future1, 3)
    result = future2.result()  # 9  (2)!
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

    # -----------------
    from ase.build import bulk

    from quacc.recipes.emt.slabs import bulk_to_defects_flow

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

    atoms = molecule("Cu")
    future = relax_job(atoms)
    result = future.result()
    assert future.done()

    # ---------------------
    from ase.build import molecule

    from quacc.recipes.lj.core import static_job

    atoms = molecule("Cu")
    future = static_job(atoms)
    result = future.result()
    assert future.done()

    # --------------------
    from ase.build import molecule

    from quacc.recipes.lj.core import freq_job

    atoms = molecule("Cu")
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
    future = static_job(
        atoms, charge=0, multiplicity=3, method="wb97m-v", basis="def2-svp"
    )
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

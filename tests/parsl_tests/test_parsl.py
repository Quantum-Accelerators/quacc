import pytest
from ase.build import bulk, molecule

try:
    import parsl
    from parsl import join_app, python_app
except ImportError:
    parsl = None


@pytest.mark.skipif(parsl is None, reason="Parsl is not installed")
def setup_module():
    parsl.load()


@pytest.mark.skipif(parsl is None, reason="Parsl is not installed")
def test_tutorial1(tmpdir):
    tmpdir.chdir()

    # Define the Python apps
    @python_app
    def relax_app(atoms):
        from quacc.recipes.emt.core import relax_job

        return relax_job(atoms)

    @python_app
    def static_app(atoms):
        from quacc.recipes.emt.core import static_job

        return static_job(atoms)

    # Make an Atoms object of a bulk Cu structure
    atoms = bulk("Cu")

    # Call App 1
    future1 = relax_app(atoms)

    # Call App 2, which takes the output of App 1 as input
    future2 = static_app(future1)
    result = future2.result()
    assert future2.done()
    assert "atoms" in result


@pytest.mark.skipif(parsl is None, reason="Parsl is not installed")
def test_tutorial2(tmpdir):
    tmpdir.chdir()

    # Define the Python app
    @python_app
    def relax_app(atoms):
        from quacc.recipes.emt.core import relax_job

        return relax_job(atoms)

    # Define two Atoms objects
    atoms1 = bulk("Cu")
    atoms2 = molecule("N2")

    # Define two independent relaxation jobs
    future1 = relax_app(atoms1)
    future2 = relax_app(atoms2)

    # Print the results
    future1.result(), future2.result()
    assert future1.done()
    assert "atoms" in future1.result()
    assert future2.done()
    assert "atoms" in future2.result()


@pytest.mark.skipif(parsl is None, reason="Parsl is not installed")
def test_tutorial3(tmpdir):
    tmpdir.chdir()

    @python_app
    def relax_app(atoms):
        from quacc.recipes.emt.core import relax_job

        return relax_job(atoms)

    @python_app
    def bulk_to_slabs_app(atoms):
        from quacc.recipes.emt.slabs import bulk_to_slabs_flow

        return bulk_to_slabs_flow(atoms, slab_static=None)

    # Define the Atoms object
    atoms = bulk("Cu")

    # Define the workflow
    future1 = relax_app(atoms)
    future2 = bulk_to_slabs_app(future1)

    # Print the results
    future2.result()
    assert len(future2.result()) == 4
    assert future2.done()


@pytest.mark.skipif(parsl is None, reason="Parsl is not installed")
def test_tutorial4(tmpdir):
    tmpdir.chdir()

    from quacc.recipes.emt.parsl.slabs import bulk_to_slabs_flow

    # Define the Python App
    @python_app
    def relax_app(atoms):
        from quacc.recipes.emt.core import relax_job

        return relax_job(atoms)

    # Define the Atoms object
    atoms = bulk("Cu")

    # Define the workflow
    future1 = relax_app(atoms)
    slab_futures = bulk_to_slabs_flow(future1, slab_static=None)

    # Print the results
    result = [slab_future.result() for slab_future in slab_futures]
    assert len(result) == 4


@pytest.mark.skipif(parsl is None, reason="Parsl is not installed")
def test_comparison1():
    @python_app
    def add(a, b):
        return a + b

    @python_app
    def mult(a, b):
        return a * b

    def workflow(a, b, c):
        return mult(add(a, b), c)

    assert workflow(1, 2, 3).result() == 9


@pytest.mark.skipif(parsl is None, reason="Parsl is not installed")
def test_comparison2():
    @python_app
    def add(a, b):
        return a + b

    @python_app
    def make_more(val):
        return [val] * 3

    @join_app
    def add_distributed(vals, c):
        return [add(val, c) for val in vals]

    def workflow(a, b, c):
        future1 = add(a, b)
        future2 = make_more(future1)
        return add_distributed(future2, c)

    assert workflow(1, 2, 3).result() == [6, 6, 6]


@pytest.mark.skipif(parsl is None, reason="Parsl is not installed")
def test_slabs(tmpdir):
    tmpdir.chdir()

    from quacc.recipes.emt.parsl.slabs import bulk_to_slabs_flow

    slab_futures = bulk_to_slabs_flow(bulk("Cu"))
    result = [slab_future.result() for slab_future in slab_futures]
    assert len(result) == 4

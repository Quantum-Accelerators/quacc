import pytest
from ase.build import bulk

try:
    from pymatgen.analysis import defects
except ImportError:
    defects = None


@pytest.mark.skipif(defects is None, reason="pymatgen-analysis-defects needed")
def test_bulk_to_defects_flow(tmpdir):
    from quacc.recipes.emt.defects import bulk_to_defects_flow

    tmpdir.chdir()

    atoms = bulk("Cu")
    output = bulk_to_defects_flow(atoms, defect_relax_kwargs={"opt_swaps": {"fmax": 5}})
    assert len(output) == 1
    assert len(output[0]["atoms"]) == 107

    atoms = bulk("Cu")
    output = bulk_to_defects_flow(
        atoms, run_static=False, defect_relax_kwargs={"opt_swaps": {"fmax": 5}}
    )
    assert len(output) == 1
    assert len(output[0]["atoms"]) == 107

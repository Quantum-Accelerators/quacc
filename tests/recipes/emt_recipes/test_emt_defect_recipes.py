import pytest

pytest.importorskip("pymatgen.analysis.defects")


def test_bulk_to_defects_flow(tmpdir):
    from ase.build import bulk

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

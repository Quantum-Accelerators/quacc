from functools import partial

import pytest
from ase.build import bulk

pytest.importorskip("pymatgen.analysis.defects")

from quacc.recipes.emt.core import relax_job
from quacc.recipes.emt.defects import bulk_to_defects_flow


def test_bulk_to_defects_flow(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    atoms = bulk("Cu")
    output = bulk_to_defects_flow(
        atoms, defect_relax_job=partial(relax_job, opt_params={"fmax": 5})
    )
    assert len(output) == 1
    assert len(output[0]["atoms"]) == 107

    atoms = bulk("Cu")
    output = bulk_to_defects_flow(
        atoms,
        defect_relax_job=partial(relax_job, opt_params={"fmax": 5}),
        defect_static_job=None,
    )
    assert len(output) == 1
    assert len(output[0]["atoms"]) == 107

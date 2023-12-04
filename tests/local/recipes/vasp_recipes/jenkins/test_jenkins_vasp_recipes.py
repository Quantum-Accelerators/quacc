from ase.build import bulk

from quacc import SETTINGS
from quacc.recipes.vasp.core import static_job

DEFAULT_SETTINGS = SETTINGS.model_copy()


def test_static_job(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    atoms = bulk("Al", kpts=[3, 3, 3])

    output = static_job(atoms)
    assert output["nsites"] == len(atoms)
    assert "isym" not in output["parameters"]
    assert output["parameters"]["nsw"] == 0
    assert output["parameters"]["lwave"] is True
    assert output["parameters"]["encut"] == 520
    assert output["parameters"]["efermi"] == "midgap"
    assert output["parameters"]["kpts"] == [3, 3, 3]

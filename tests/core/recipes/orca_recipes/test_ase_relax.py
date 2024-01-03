import os
from shutil import which

import pytest
from ase.build import molecule

from quacc import SETTINGS
from quacc.recipes.orca.core import ase_relax_job

orca_path = which(SETTINGS.ORCA_CMD)

has_orca = bool(orca_path and os.path.getsize(orca_path) > 1024 * 1024)
pytestmark = pytest.mark.skipif(not has_orca, reason="Needs ORCA")


def test_ase_relax_job(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    atoms = molecule("H2")

    output = ase_relax_job(atoms, opt_params={"fmax": 0.1}, nprocs=1)
    assert output["natoms"] == len(atoms)
    assert output["parameters"]["charge"] == 0
    assert output["parameters"]["mult"] == 1
    assert (
        output["parameters"]["orcasimpleinput"]
        == "wb97x-d3bj def2-tzvp engrad slowconv normalprint xyzfile"
    )
    assert output.get("trajectory") is not None
    assert output.get("trajectory_results")
    assert output.get("attributes") is not None

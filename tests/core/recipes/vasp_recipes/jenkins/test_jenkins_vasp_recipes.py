from __future__ import annotations

from shutil import which

import pytest

from quacc import get_settings

pytestmark = pytest.mark.skipif(
    which(str(get_settings().VASP_CMD)) is None, reason="VASP not installed"
)

import numpy as np
from ase.build import bulk
from numpy.testing import assert_equal
from pymatgen.core import Structure

from quacc.recipes.vasp.core import static_job


def test_static_job(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    atoms = bulk("Al")

    output = static_job(atoms, kpts=[3, 3, 3])
    assert output["structure_metadata"]["nsites"] == len(atoms)
    assert "isym" not in output["parameters"]
    assert output["parameters"]["nsw"] == 0
    assert output["parameters"]["lwave"] is True
    assert output["parameters"]["encut"] == 520
    assert output["parameters"]["efermi"] == "midgap"
    assert output["parameters"]["kpts"] == [3, 3, 3]
    assert output["parameters"]["pp"] == "PBE"
    assert output["parameters"]["pp_version"] == "64"
    assert output["results"]["energy"] < 0


def test_static_job_spin(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    atoms = bulk("Fe")
    atoms.set_initial_magnetic_moments([5.0] * len(atoms))

    output = static_job(atoms, kpts=[3, 3, 3])
    assert output["structure_metadata"]["nsites"] == len(atoms)
    assert "isym" not in output["parameters"]
    assert output["parameters"]["nsw"] == 0
    assert output["parameters"]["lwave"] is True
    assert output["parameters"]["encut"] == 520
    assert output["parameters"]["efermi"] == "midgap"
    assert output["parameters"]["kpts"] == [3, 3, 3]
    assert output["results"]["energy"] < 0
    output_magmoms = np.array(
        Structure.from_dict(output["structure"]).site_properties["magmom"]
    )
    assert output_magmoms.all()
    assert_equal(output["atoms"].get_initial_magnetic_moments(), output_magmoms)

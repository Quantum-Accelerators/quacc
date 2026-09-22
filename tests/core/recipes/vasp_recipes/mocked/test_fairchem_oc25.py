from __future__ import annotations

import importlib.util
from copy import deepcopy
from pathlib import Path

import numpy as np
import pytest
from ase import Atoms
from ase.constraints import FixAtoms
from monty.io import zopen
from monty.os.path import zpath
from pymatgen.io.vasp import Incar, Kpoints

from quacc.calculators.vasp import Vasp
from quacc.recipes.vasp import fairchem

requires_fairchem_oc = pytest.mark.skipif(
    not fairchem.has_fairchem_oc, reason="fairchem-data-oc not installed"
)


@pytest.fixture
def interface(tmp_path, monkeypatch):
    pseudo_dir = tmp_path / "pseudos"
    for setup in ("Na", "Na_pv", "O", "H"):
        potential = pseudo_dir / "potpaw_PBE.64" / setup / "POTCAR"
        potential.parent.mkdir(parents=True)
        potential.write_text("<3\n")
    monkeypatch.setenv("VASP_PP_PATH", str(pseudo_dir))
    atoms = Atoms(
        "NaOH",
        positions=[[0, 0, 8], [0, 0, 10], [0.9, 0, 10.4]],
        cell=[8, 8, 24],
        pbc=True,
        tags=[0, 2, 2],
        constraint=FixAtoms(indices=[0]),
    )
    atoms.info["sid"] = "oc25-test"
    return atoms


@requires_fairchem_oc
@pytest.mark.parametrize("skew", [False, True])
def test_fairchem_oc25_inputs(
    patch_nonmetallic_taskdoc, interface, tmp_path, monkeypatch, skew
):
    from fairchem.data.oc.utils.vasp import write_vasp_input_files
    from fairchem.data.oc.utils.vasp_flags import SOLVENT_BASE_FLAGS

    if skew:
        interface.set_cell([[8, 6, 0], [0, 8, 0], [0, 0, 24]])
    original = interface.copy()
    original_flags = deepcopy(SOLVENT_BASE_FLAGS)
    observed = {}
    write_input = Vasp.write_input

    def capture_input(self, *args, **kwargs):
        write_input(self, *args, **kwargs)
        observed["potentials"] = [Path(p) for p in self.ppp_list]
        observed["copilot"] = self.incar_copilot_mode
        observed["custodian"] = self.use_custodian

    monkeypatch.setattr(Vasp, "write_input", capture_input)
    output = fairchem.oc25_static_job(interface)
    directory = Path(output["dir_name"])

    reference = tmp_path / "reference"
    reference.mkdir()
    write_vasp_input_files(
        interface.copy(),
        reference,
        vasp_flags=SOLVENT_BASE_FLAGS
        | {"ediff": 1e-6, "nsw": 0, "ibrion": -1, "pp_version": "64"},
        pp_setups="recommended",
    )

    actual_incar = Incar.from_file(zpath(str(directory / "INCAR")))
    reference_incar = Incar.from_file(reference / "INCAR")
    assert actual_incar == reference_incar
    assert actual_incar["GGA"].upper() == "RP"
    assert actual_incar["IVDW"] == 11
    assert actual_incar["NSW"] == 0
    assert actual_incar["IBRION"] == -1
    assert actual_incar["EDIFF"] == 1e-6
    assert actual_incar["ISPIN"] == 1

    actual_kpoints = Kpoints.from_file(zpath(str(directory / "KPOINTS")))
    assert actual_kpoints == Kpoints.from_file(reference / "KPOINTS")
    assert actual_kpoints.kpts == [(5, 5, 1)]
    assert actual_kpoints.style == Kpoints.supported_modes.Monkhorst
    for filename in ("POSCAR", "POTCAR"):
        with zopen(zpath(str(directory / filename)), mode="rt", encoding="utf-8") as f:
            assert f.read() == (reference / filename).read_text()

    # The fake POTCAR contents are identical: check the actual selected paths too.
    assert [p.parent.name for p in observed["potentials"]] == ["Na_pv", "O", "H"]
    assert all(p.parent.parent.name == "potpaw_PBE.64" for p in observed["potentials"])
    assert observed["copilot"] == "off"
    assert observed["custodian"] is False
    assert output["name"] == "OC25 Static"
    assert output["structure_metadata"]["nsites"] == len(interface)
    assert np.shape(output["results"]["forces"]) == (len(interface), 3)
    assert np.isfinite(output["results"]["energy"])
    assert interface == original
    assert interface.info == original.info
    np.testing.assert_array_equal(interface.get_tags(), original.get_tags())
    assert interface.constraints[0].todict() == original.constraints[0].todict()
    assert original_flags == SOLVENT_BASE_FLAGS


@requires_fairchem_oc
def test_fairchem_oc25_overrides(patch_nonmetallic_taskdoc, interface, tmp_path):
    (tmp_path / "seed.txt").write_text("copied input")
    output = fairchem.oc25_static_job(
        interface,
        ediff=1e-4,
        encut=450,
        kpts=(3, 2, 1),
        ncore=2,
        dipol=None,
        additional_fields={"name": "Custom OC25", "source_id": "test-frame"},
        copy_files=[{"source": tmp_path, "filenames": ["seed.txt"]}],
    )
    assert output["parameters"]["ediff"] == 1e-4
    assert output["parameters"]["encut"] == 450
    assert output["parameters"]["kpts"] == (3, 2, 1)
    assert output["parameters"]["ncore"] == 2
    assert "dipol" not in output["parameters"]
    assert output["parameters"]["nsw"] == 0
    assert output["name"] == "Custom OC25"
    assert output["source_id"] == "test-frame"
    with zopen(
        zpath(str(Path(output["dir_name"]) / "seed.txt")), mode="rt", encoding="utf-8"
    ) as f:
        assert f.read() == "copied input"


def test_fairchem_oc25_missing_dependency(monkeypatch, interface):
    find_spec = importlib.util.find_spec

    def without_fairchem(name, *args, **kwargs):
        return None if name == "fairchem" else find_spec(name, *args, **kwargs)

    # Evaluate the normal dependency guard without reloading the shared module.
    spec = importlib.util.spec_from_file_location(
        "fairchem_without_oc", fairchem.__file__
    )
    module = importlib.util.module_from_spec(spec)
    monkeypatch.setattr(importlib.util, "find_spec", without_fairchem)
    spec.loader.exec_module(module)
    with pytest.raises(RuntimeError, match=r"fairchem-data-oc.*quacc\[fairchem\]"):
        module.oc25_static_job(interface)

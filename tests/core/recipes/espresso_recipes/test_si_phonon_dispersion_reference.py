"""Regression checks against unmodified QE 7.5 matdyn.x output for Si.

The reference files were produced from the endpoint list in the manifest,
passed to matdyn.x with ``q_in_band_form=True`` and ``q_in_cryst_coord=True``.

TODO: Si has no non-analytic term. A small polar reference (e.g. NaCl or h-BN)
would cover the LO-TO limits at Gamma in tests that require no QE executable.
"""

from __future__ import annotations

import gzip
import hashlib
import json
from io import StringIO
from pathlib import Path

import numpy as np
import pytest
from ase import Atoms
from ase.units import Bohr
from numpy.testing import assert_allclose, assert_array_equal

from quacc.calculators.espresso.phonons import (
    _check_matdyn_qpoints,
    prepare_phonon_path,
)
from quacc.calculators.espresso.utils import read_matdyn_freq

DATA = Path(__file__).parent / "data" / "si_phonon_dispersion"


def reference():
    manifest = json.loads((DATA / "manifest.json").read_text())
    files = {}
    for name, digest in manifest["output_sha256"].items():
        content = gzip.decompress((DATA / (name + ".gz")).read_bytes())
        assert hashlib.sha256(content).hexdigest() == digest
        files[name] = content.decode()
    return manifest, files


@pytest.mark.parametrize("layout", ["freq", "gp", "both"])
def test_si_native_frequency_files(tmp_path, layout):
    manifest, files = reference()
    if layout != "gp":
        (tmp_path / "matdyn.freq").write_text(files["matdyn.freq"])
    if layout != "freq":
        (tmp_path / "matdyn.freq.gp").write_text(files["matdyn.freq.gp"])
    parsed = read_matdyn_freq(tmp_path / "matdyn.freq")
    gp = np.loadtxt(StringIO(files["matdyn.freq.gp"]))
    # With six Si modes, each native frequency record occupies one line.
    native = np.array(
        [
            np.fromstring(line, sep=" ")
            for line in files["matdyn.freq"].splitlines()[2::2]
        ]
    )
    assert_array_equal(parsed["frequencies"], native)
    assert_array_equal(parsed["frequencies"], gp[:, 1:])
    assert_array_equal(np.signbit(parsed["frequencies"]), np.signbit(native))
    assert parsed["frequencies"].shape == tuple(manifest["frequency_shape"])
    assert parsed["frequency_units"] == "cm^-1"
    if layout == "gp":
        assert "qpoints" not in parsed
    else:
        assert parsed["qpoints"].shape == (len(gp), 3)
        assert parsed["qpoint_units"] == "2pi/alat"
    if layout == "freq":
        assert "distances" not in parsed
    else:
        assert_array_equal(parsed["distances"], gp[:, 0])


def test_si_endpoints_ticks_and_basis(tmp_path):
    manifest, files = reference()
    atoms = Atoms("Si2", cell=manifest["cell_angstrom"], pbc=True)
    endpoints = np.array(manifest["matdyn_qpts"])
    path = prepare_phonon_path(
        atoms, False, 20, endpoints, manifest["tick_labels"], None
    )
    npoints = 1 + int(np.maximum(endpoints[:-1, 3], 1).sum())
    assert len(path["qpoints"]) == npoints == manifest["frequency_shape"][0]
    assert_allclose(
        path["qpoints"][path["tick_indices"]], endpoints[:, :3], rtol=0, atol=1e-14
    )
    assert path["tick_indices"] == manifest["tick_indices"]
    assert [list(segment) for segment in path["segments"]] == manifest["segments"]
    assert_allclose(
        path["tick_positions"], manifest["tick_positions"], rtol=0, atol=1e-12
    )
    for name, content in files.items():
        (tmp_path / name).write_text(content)
    parsed = read_matdyn_freq(tmp_path / "matdyn.freq")
    alat_angstrom = manifest["alat_bohr"] * Bohr
    cartesian = path["qpoints"] @ atoms.cell.reciprocal().array * alat_angstrom
    assert_allclose(parsed["qpoints"], cartesian, rtol=0, atol=2e-6)
    distances = parsed["distances"].copy()
    # Remove the jumps recorded by QE's .gp writer at the native zero counts.
    for start, _stop in path["segments"][1:]:
        distances[start:] -= distances[start] - distances[start - 1]
    assert_allclose(
        path["distances"], distances * 2 * np.pi / alat_angstrom, rtol=0, atol=5e-6
    )
    assert_allclose(
        path["tick_positions"], path["distances"][path["tick_indices"]], rtol=0, atol=0
    )


def test_si_basis_check_with_qe_output(tmp_path):
    """matdyn.x's Cartesian q-points identify the cell used for the path."""
    manifest, files = reference()
    (tmp_path / "matdyn.freq").write_text(files["matdyn.freq"])
    cartesian = read_matdyn_freq(tmp_path / "matdyn.freq")["qpoints"]
    cell = np.array(manifest["cell_angstrom"])
    path = prepare_phonon_path(
        Atoms("Si2", cell=cell, pbc=True),
        False,
        20,
        manifest["matdyn_qpts"],
        manifest["tick_labels"],
        None,
    )
    # Gamma starts the second segment, after a break from another point.
    assert path["segments"][1][0] == 131
    assert path["tick_labels"][path["tick_indices"].index(131)] == "Γ"
    assert path["matdyn_indices"][131] == 132
    # The reference was run without the auxiliary point, so compare the path itself.
    _check_matdyn_qpoints({"matdyn_qpts": path["qpoints"], "cell": cell}, cartesian)
    with pytest.raises(ValueError, match="inconsistent with the cell"):
        _check_matdyn_qpoints(
            {"matdyn_qpts": path["qpoints"], "cell": cell[[1, 0, 2]]}, cartesian
        )


def test_si_sampling_inputs():
    manifest, _ = reference()
    ph = manifest["inputs"]["ph"]
    assert [ph[key] for key in ("nq1", "nq2", "nq3")] == [3, 3, 3]
    assert ph["ldisp"] is True
    assert ph["qplot"] is False
    assert manifest["ph_cards"] == []
    matdyn = manifest["inputs"]["matdyn"]
    assert matdyn["dos"] is False
    assert matdyn["q_in_band_form"] is True
    assert matdyn["q_in_cryst_coord"] is True
    assert not {"nk1", "nk2", "nk3"}.intersection(matdyn)

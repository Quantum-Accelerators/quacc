from __future__ import annotations

from io import StringIO
from typing import TYPE_CHECKING
from unittest.mock import Mock

import numpy as np
from ase import Atoms
from ase.io.espresso import write_espresso_in
from numpy.testing import assert_allclose

from quacc.recipes.espresso.bands import bands_pw_job
from quacc.utils.kpts import convert_pmg_kpts

if TYPE_CHECKING:
    import pytest


def test_bands_pw_job_orthorhombic_kpoints(monkeypatch: pytest.MonkeyPatch) -> None:
    atoms = Atoms("Al", cell=[5.0, 4.0, 3.2], pbc=True)
    expected_kpts = convert_pmg_kpts({"line_density": 1}, atoms, force_gamma=True)[0]
    run_and_summarize = Mock()
    monkeypatch.setattr(
        "quacc.recipes.espresso.bands.run_and_summarize", run_and_summarize
    )

    bands_pw_job(atoms, line_density=1)

    assert run_and_summarize.call_args.args[0] is atoms
    path = run_and_summarize.call_args.kwargs["calc_defaults"]["kpts"]
    pw_input = StringIO()
    write_espresso_in(pw_input, atoms, kpts=path, pseudopotentials={"Al": "Al.upf"})
    lines = pw_input.getvalue().splitlines()
    start = lines.index("K_POINTS crystal_b") + 1
    assert int(lines[start]) == len(expected_kpts)
    written_kpts = np.loadtxt(lines[start + 1 : start + 1 + len(expected_kpts)])

    assert_allclose(path.cell, atoms.get_cell())
    assert_allclose(path.kpts, expected_kpts, atol=1.0e-14, rtol=0)
    assert_allclose(written_kpts[:, :3], expected_kpts, atol=1.0e-14, rtol=0)
    assert_allclose(written_kpts[:, 3], 0)
    assert_allclose(
        path.cartesian_kpts(),
        expected_kpts @ atoms.cell.reciprocal(),
        atol=1.0e-14,
        rtol=0,
    )

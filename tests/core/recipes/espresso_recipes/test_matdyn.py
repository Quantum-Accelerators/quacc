from __future__ import annotations

from pathlib import Path
from shutil import which

import numpy as np
import pytest
from ase.io.espresso import read_fortran_namelist
from monty.io import zopen
from numpy.testing import assert_allclose, assert_array_equal

from quacc import change_settings, get_settings
from quacc.calculators.espresso.espresso import EspressoTemplate
from quacc.recipes.espresso.phonons import matdyn_job

DATA_DIR = Path(__file__).parent / "data" / "matdyn_test"
MATDYN = str(
    get_settings().ESPRESSO_BIN_DIR / get_settings().ESPRESSO_BINARIES["matdyn"]
)


def test_matdyn_job_structured_results(tmp_path, monkeypatch):
    """Exercise the job, staging, input writer, parser, and schema without QE."""
    monkeypatch.chdir(tmp_path)
    source = tmp_path / "source"
    source.mkdir()
    (source / "q2r.fc").write_text("mock force constants")

    def execute(self, directory, profile):
        assert (directory / "q2r.fc").read_text() == "mock force constants"
        with (directory / "matdyn.in").open() as handle:
            inputs, cards = read_fortran_namelist(handle)
        assert inputs["input"]["q_in_cryst_coord"] is True
        assert inputs["input"]["q_in_band_form"] is False
        assert_array_equal(np.loadtxt(cards[1:-1]), [[0, 0, 0], [0.5, 0, 0]])
        assert inputs["input"]["flfrq"] == "custom.freq"
        (directory / "custom.freq").write_text(
            " &plot nbnd=3, nks=2 /\n0 0 0\n-1 2 3\n0.5 0.5 0\n4 5 6\n"
        )
        (directory / "custom.freq.gp").write_text("0 -1 2 3\n0.707107 4 5 6\n")

    monkeypatch.setattr(EspressoTemplate, "execute", execute)
    result = matdyn_job(
        [{"source": source, "filenames": "q2r.fc"}],
        {"purpose": "dispersion"},  # Keep the existing positional API.
        qpts=[[0, 0, 0], [0.5, 0, 0]],
        input_data={
            "input": {
                "dos": False,
                "q_in_band_form": False,
                "q_in_cryst_coord": True,
                "flfrq": "custom.freq",
            }
        },
    )
    assert result["purpose"] == "dispersion"
    assert result["name"] == "matdyn.x Phonon"
    assert_array_equal(result["parameters"]["qpts"], [[0, 0, 0], [0.5, 0, 0]])
    parsed = result["results"]["matdyn_results"]
    assert_array_equal(parsed["frequencies"], [[-1, 2, 3], [4, 5, 6]])
    assert parsed["frequency_units"] == "cm^-1"
    assert parsed["qpoint_coordinate_system"] == "cartesian"
    assert_array_equal(parsed["qpoints"], [[0, 0, 0], [0.5, 0.5, 0]])
    assert_array_equal(parsed["distances"], [0, 0.707107])


@pytest.mark.skipif(which(MATDYN) is None, reason="matdyn.x not installed")
@pytest.mark.parametrize("path_type", ["native", "band", "explicit", "disconnected"])
def test_matdyn_dispersion_qe(tmp_path, monkeypatch, path_type):
    """Check QE's actual interpolation counts with the existing IFC fixture.

    This tests executable interoperability, not convergence or physical accuracy.
    """
    monkeypatch.chdir(tmp_path)
    monkeypatch.setenv("OMP_NUM_THREADS", "1")
    points = np.array(
        [
            [0.5, 0.5, 0.5, 20],
            [0, 0, 0, 20],
            [0.5, 0, 0.5, 20],
            [0.375, 0.375, 0.75, 20],
            [0, 0, 0, 1],
        ]
    )
    indices = [0, 20, 40, 60, 80]
    nq = 81
    if path_type == "explicit":
        points = points[:, :3]
        indices = list(range(5))
        nq = 5
    elif path_type == "disconnected":
        points = np.array(
            [[0, 0, 0, 2], [0.5, 0, 0, 0], [0.5, 0.5, 0, 3], [0, 0, 0, 1]]
        )
        indices = [0, 2, 3, 6]
        nq = 7
    kwargs = {"qpts": points}
    if path_type == "native":
        kwargs = {
            "additional_cards": ["5"]
            + [" ".join(map(str, point[:3])) + f" {int(point[3])}" for point in points]
        }

    with change_settings({"ESPRESSO_PARALLEL_CMD": ""}):
        result = matdyn_job(
            [{"source": DATA_DIR, "filenames": "q2r.fc"}],
            input_data={
                "input": {
                    "dos": False,
                    "q_in_band_form": path_type != "explicit",
                    "q_in_cryst_coord": True,
                }
            },
            **kwargs,
        )
    parsed = result["results"]["matdyn_results"]
    assert parsed["frequencies"].shape == (nq, 3)
    assert parsed["qpoints"].shape == (nq, 3)
    # Reciprocal basis of the bcc lattice in q2r.fc, expressed in QE's units.
    direct = np.array([[-1, 1, 1], [1, -1, 1], [1, 1, -1]]) / np.sqrt(3)
    expected = points[:, :3] @ np.linalg.inv(direct).T
    assert_allclose(parsed["qpoints"][indices], expected, atol=1e-6)
    with zopen(Path(result["dir_name"]) / "matdyn.in.gz", "rt") as handle:
        inputs = read_fortran_namelist(handle)[0]["input"]
    assert not {"nk1", "nk2", "nk3"}.intersection(inputs)

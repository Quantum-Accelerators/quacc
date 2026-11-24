from __future__ import annotations

import gzip
from copy import deepcopy
from pathlib import Path
from shutil import which

import numpy as np
import pytest
from ase import Atoms
from ase.io.espresso import read_fortran_namelist
from numpy.testing import assert_allclose, assert_array_equal

from quacc import change_settings, get_settings
from quacc.calculators.espresso.phonons import (
    prepare_phonon_path,
    summarize_phonon_dispersion,
)
from quacc.recipes.espresso._base import prepare_calc
from quacc.recipes.espresso.core import static_job
from quacc.recipes.espresso.phonons import _PHONON_ASR_DEFAULTS, phonon_dispersion_flow
from quacc.utils.files import copy_decompress_files


class StagedCalls(dict):
    """Calls made to each mocked binary.

    ``qe_cell`` is the cell of the mocked force constants, which matdyn.x uses
    to convert fractional q-points to the Cartesian ones it writes. ``has_zstar``
    sets the dielectric flag written by the mocked q2r.x.
    """

    qe_cell = np.diag([2.0, 3.0, 4.0])
    has_zstar = True


def force_constants(has_zstar):
    """Header of a q2r.x text force-constant file for one H atom (ibrav=0)."""
    return (
        "  1    1   0  2.0000000  0.0000000  0.0000000  0.0000000\n"
        "      1.000000000    0.000000000    0.000000000\n"
        "      0.000000000    1.500000000    0.000000000\n"
        "      0.000000000    0.000000000    2.000000000\n"
        "           1  'H     '    918.68\n"
        "    1    1      0.0000000000      0.0000000000      0.0000000000\n"
        f" {'T' if has_zstar else 'F'}   1.2\n"
    )


@pytest.fixture(params=[False, True], ids=["text", "gzip"])
def staged_dispersion(tmp_path, monkeypatch, request):
    """Mock executable output while using the real input writer, parser, and transfers."""
    calls = StagedCalls()

    def write_file(directory, name, content):
        if request.param:
            with gzip.open(directory / (name + ".gz"), "wt") as handle:
                handle.write(content)
        else:
            (directory / name).write_text(content)

    def run_and_summarize(
        *, template, calc_defaults, calc_swaps, copy_files, additional_fields
    ):
        binary = template.binary
        directory = tmp_path / binary
        directory.mkdir()
        for spec in copy_files or []:
            copy_decompress_files(spec["source"], spec["filenames"], directory)
        calc = prepare_calc(
            Atoms(),
            template=template,
            calc_defaults=calc_defaults,
            calc_swaps=calc_swaps,
        )
        calc.directory = directory
        calc.write_inputfiles(Atoms(), ["energy"])
        calls[binary] = {
            "parameters": deepcopy(calc.parameters),
            "copy_files": deepcopy(copy_files),
            "prev_outdir": template.outdir,
        }
        if binary == "ph":
            write_file(directory, "matdyn0", "q grid")
            write_file(directory, "matdyn1", "dynamical matrix")
            # XML dynamical matrices and restart files must stay distinguishable.
            write_file(directory, "matdyn2.xml", "xml dynamical matrix")
            (directory / "pwscf.save").mkdir()
            (directory / "pwscf.save/wfc1.dat").write_text("wavefunctions")
        elif binary == "q2r":
            assert (directory / "matdyn0").read_text() == "q grid"
            assert (directory / "matdyn1").read_text() == "dynamical matrix"
            assert (directory / "matdyn2.xml").read_text() == "xml dynamical matrix"
            write_file(directory, "q2r.fc", force_constants(calls.has_zstar))
            write_file(directory, "q2r.fc.xml", "xml force constants")
        else:
            fc_text = (directory / "q2r.fc").read_text()
            assert fc_text == force_constants(calls.has_zstar)
            assert (directory / "q2r.fc.xml").read_text() == "xml force constants"
            with (directory / "matdyn.in").open() as handle:
                inputs, cards = read_fortran_namelist(handle)
            assert not inputs["input"].get("q_in_band_form", False)
            points = np.loadtxt(cards[1:-1], ndmin=2)
            assert points.shape == (int(cards[0]), 3)
            # matdyn.x writes Cartesian q-points in units of 2pi/alat.
            alat = np.linalg.norm(calls.qe_cell[0])
            cartesian = points @ np.linalg.inv(calls.qe_cell).T * alat
            frequencies = np.arange(3.0 * len(points)).reshape(-1, 3)
            frequencies[0, 0] = -5
            content = f" &plot nbnd=3, nks={len(points)} /\n"
            for point, row in zip(cartesian, frequencies, strict=True):
                content += " ".join(f"{value:.6f}" for value in point) + "\n"
                content += " ".join(map(str, row)) + "\n"
            filename = inputs["input"]["flfrq"]
            (directory / filename).write_text(content)
        results = template.read_results(directory) if binary == "matdyn" else {}
        write_file(directory, f"{binary}.out", "mock executable output")
        return {
            "dir_name": directory,
            "parameters": calc.parameters,
            "results": results,
            **additional_fields,
        }

    monkeypatch.setattr(
        "quacc.recipes.espresso.phonons.run_and_summarize", run_and_summarize
    )
    return calls


@pytest.mark.parametrize("disconnected", [False, True])
def test_dispersion_flow_staging(tmp_path, staged_dispersion, disconnected):
    scf = tmp_path / "scf"
    scf.mkdir()
    (scf / "scf.dat").write_text("SCF input")
    atoms = Atoms("H", cell=[2, 3, 4], pbc=True)
    points = [
        [0, 0, 0, 2],
        [0.5, 0, 0, 0 if disconnected else 2],
        [0, 0.5, 0, 3],
        [0, 0.5, 0.5, 1],
    ]
    params = {
        "all": {"additional_fields": {"tag": "dispersion"}},
        "phonon_job": {
            "input_data": {"inputph": {"nq1": 2, "nq2": 3, "nq3": 4, "epsil": True}}
        },
        "q2r_job": {"input_data": {"input": {"write_lr": True}}},
        "matdyn_job": {
            "input_data": {
                "input": {
                    "asr": "all",
                    "read_lr": True,
                    "loto_disable": False,
                    "flfrq": "custom.freq",
                }
            }
        },
    }
    original = deepcopy(params)
    defaults = deepcopy(_PHONON_ASR_DEFAULTS)
    copy_files = [{"source": scf, "filenames": "scf.dat"}]
    result = phonon_dispersion_flow(
        atoms,
        copy_files=copy_files,
        prev_outdir=scf,
        qpts=points,
        qpoint_labels=["G", "X", "Y", "Z"],
        job_params=params,
    )
    assert set(result) == {"phonon_job", "q2r_job", "matdyn_job", "dispersion"}
    assert list(staged_dispersion) == ["ph", "q2r", "matdyn"]
    assert staged_dispersion["ph"]["copy_files"] == copy_files
    assert staged_dispersion["ph"]["prev_outdir"] == scf
    ph_input = staged_dispersion["ph"]["parameters"]["input_data"]["inputph"]
    assert ph_input["ldisp"] is True
    assert ph_input["qplot"] is False
    assert [ph_input[key] for key in ("nq1", "nq2", "nq3")] == [2, 3, 4]
    assert ph_input["epsil"] is True
    matdyn_input = staged_dispersion["matdyn"]["parameters"]["input_data"]["input"]
    assert matdyn_input["dos"] is False
    assert matdyn_input["q_in_cryst_coord"] is True
    assert not matdyn_input.get("q_in_band_form", False)
    assert matdyn_input["loto_disable"] is False
    assert matdyn_input["read_lr"] is True
    assert not {"nk1", "nk2", "nk3"}.intersection(matdyn_input)
    assert not (tmp_path / "q2r/scf.dat").exists()
    assert not (tmp_path / "q2r/pwscf.save").exists()
    assert not list((tmp_path / "matdyn").glob("matdyn[0-9]*"))
    assert not list((tmp_path / "matdyn").glob("q2r.out*"))
    assert (tmp_path / "ph/pwscf.save/wfc1.dat").exists()
    for name in ("phonon_job", "q2r_job", "matdyn_job"):
        assert result[name]["tag"] == "dispersion"
    dispersion = result["dispersion"]
    nq = 7 if disconnected else 8
    assert dispersion["qpoints"].shape == (nq, 3)
    assert dispersion["frequencies"].shape == (nq, 3)
    assert dispersion["frequencies"][0, 0] == -5
    assert dispersion["frequency_units"] == "cm^-1"
    assert dispersion["qpoint_coordinate_system"] == "fractional"
    assert dispersion["distance_units"] == "angstrom^-1"
    assert dispersion["segments"] == ([(0, 3), (3, 7)] if disconnected else [(0, 8)])
    if disconnected:
        assert dispersion["distances"][2] == dispersion["distances"][3]
    assert dispersion["tick_labels"] == ["G", "X", "Y", "Z"]
    assert_array_equal(
        dispersion["tick_positions"],
        dispersion["distances"][dispersion["tick_indices"]],
    )
    assert dispersion["metadata"]["qgrid"] == [2, 3, 4]
    assert dispersion["metadata"]["asr"] == "all"
    assert dispersion["metadata"]["zasr"] == "crystal"
    assert dispersion["metadata"]["loto_disable"] is False
    assert dispersion["metadata"]["loto_2d"] is False
    assert dispersion["metadata"]["has_zstar"] is True
    assert_allclose(dispersion["metadata"]["cell"], atoms.cell)
    assert params == original
    assert defaults == _PHONON_ASR_DEFAULTS


def test_gamma_segment_start_frequencies_discarded(tmp_path, staged_dispersion):
    """The auxiliary Gamma evaluated before a segment is not returned."""
    result = phonon_dispersion_flow(
        Atoms("H", cell=[2, 3, 4], pbc=True),
        prev_outdir=tmp_path,
        qpts=[[0.5, 0, 0, 1], [0, 0.5, 0, 0], [0, 0, 0, 2], [0, 0, 0.5, 1]],
        qpoint_labels=["X", "Y", "G", "Z"],
    )
    raw = result["matdyn_job"]["results"]["matdyn_results"]["frequencies"]
    assert raw.shape == (6, 3)
    dispersion = result["dispersion"]
    assert dispersion["segments"] == [(0, 2), (2, 5)]
    assert_array_equal(dispersion["frequencies"], raw[[0, 1, 3, 4, 5]])
    assert_array_equal(dispersion["qpoints"][2], [0, 0, 0])
    assert dispersion["tick_labels"] == ["X", "Y", "G", "Z"]


def test_manual_explicit_and_decorators(tmp_path, staged_dispersion):
    decorated = []

    def decorator(func):
        def wrapper(*args, **kwargs):
            decorated.append(func.__name__)
            return func(*args, **kwargs)

        return wrapper

    result = phonon_dispersion_flow(
        Atoms("H", cell=[2, 3, 4], pbc=True),
        prev_outdir=tmp_path,
        make_bandpath=False,
        qpts=[[0, 0, 0], [0.5, 0, 0], [0, 0.5, 0], [0, 0.5, 0.5]],
        qpoint_labels=["G", "X", "Y", "Z"],
        segment_starts=[0, 2],
        job_decorators={"all": decorator},
    )
    assert decorated == ["phonon_job", "q2r_job", "matdyn_job"]
    assert result["dispersion"]["segments"] == [(0, 2), (2, 4)]
    assert result["dispersion"]["frequencies"].shape == (4, 3)
    assert result["dispersion"]["metadata"]["asr"] == "crystal"
    assert result["dispersion"]["metadata"]["qgrid"] == [4, 4, 4]


def test_automatic_dispersion_flow(tmp_path, staged_dispersion):
    atoms = Atoms("H", cell=[5, 4, 3.2], pbc=True)
    staged_dispersion.qe_cell = atoms.cell.array
    result = phonon_dispersion_flow(atoms, prev_outdir=tmp_path, line_density=1)
    dispersion = result["dispersion"]
    path = prepare_phonon_path(atoms, True, 1, None, None, None)
    raw = result["matdyn_job"]["results"]["matdyn_results"]["frequencies"]
    assert dispersion["metadata"]["path_source"] == "latimer_munro"
    assert dispersion["tick_labels"] == path["tick_labels"]
    assert dispersion["segments"] == path["segments"]
    assert_array_equal(dispersion["frequencies"], raw[path["matdyn_indices"]])
    assert dispersion["qpoints"].shape == dispersion["frequencies"].shape
    assert dispersion["distances"].shape == (len(dispersion["qpoints"]),)


def test_scf_cell_mismatch(tmp_path, staged_dispersion):
    """A path built for a differently oriented cell must not be mislabeled."""
    staged_dispersion.qe_cell = np.diag([3.0, 2.0, 4.0])
    with pytest.raises(ValueError, match="inconsistent with the cell of atoms"):
        phonon_dispersion_flow(
            Atoms("H", cell=[2, 3, 4], pbc=True),
            prev_outdir=tmp_path,
            qpts=[[0, 0, 0, 2], [0.5, 0, 0, 2], [0.5, 0.5, 0, 1]],
        )


def test_alat_independent_basis_check(tmp_path, staged_dispersion):
    """Only the reciprocal basis is compared; QE's alat convention may differ."""
    staged_dispersion.qe_cell = 1.01 * np.diag([2.0, 3.0, 4.0])
    result = phonon_dispersion_flow(
        Atoms("H", cell=[2, 3, 4], pbc=True),
        prev_outdir=tmp_path,
        qpts=[[0, 0, 0, 2], [0.5, 0, 0, 2], [0.5, 0.5, 0, 1]],
    )
    assert result["dispersion"]["frequencies"].shape == (5, 3)


def test_numpy_qgrid(tmp_path, staged_dispersion):
    qgrid = {f"nq{i}": np.int64(2) for i in (1, 2, 3)}
    result = phonon_dispersion_flow(
        Atoms("H", cell=[2, 3, 4], pbc=True),
        prev_outdir=tmp_path,
        qpts=[[0, 0, 0], [0.5, 0, 0]],
        job_params={"phonon_job": {"input_data": {"inputph": qgrid}}},
    )
    assert result["dispersion"]["metadata"]["qgrid"] == [2, 2, 2]


@pytest.mark.parametrize(
    ("has_zstar", "inputs"),
    [
        (True, {}),
        (False, {}),
        (True, {"loto_disable": True}),
        (True, {"loto_2d": True}),
    ],
)
def test_lo_to_metadata(tmp_path, staged_dispersion, has_zstar, inputs):
    """The dispersion records the inputs and data that control LO-TO splitting."""
    staged_dispersion.has_zstar = has_zstar
    params = {"matdyn_job": {"input_data": {"input": inputs}}}
    if "loto_2d" in inputs:
        params["q2r_job"] = {"input_data": {"input": {"loto_2d": True}}}
    result = phonon_dispersion_flow(
        Atoms("H", cell=[2, 3, 4], pbc=True),
        prev_outdir=tmp_path,
        qpts=[[0.5, 0, 0], [0, 0, 0]],
        job_params=params,
    )
    metadata = result["dispersion"]["metadata"]
    assert metadata["has_zstar"] is has_zstar
    assert metadata["loto_disable"] is inputs.get("loto_disable", False)
    assert metadata["loto_2d"] is inputs.get("loto_2d", False)
    parsed = result["matdyn_job"]["results"]["matdyn_results"]
    assert parsed["has_zstar"] is has_zstar


@pytest.mark.parametrize(
    ("params", "message"),
    [
        ({"phonon_job": {"input_data": {"ldisp": False}}}, "ldisp=True"),
        ({"phonon_job": {"input_data": {"qplot": True}}}, "qplot=False"),
        ({"phonon_job": {"use_phcg": True}}, "ph.x"),
        ({"phonon_job": {"input_data": {"nq1": 0}}}, "positive integers"),
        ({"phonon_job": {"input_data": {"nq1": 2.0}}}, "positive integers"),
        ({"phonon_job": {"input_data": {"nq1": True}}}, "positive integers"),
        ({"matdyn_job": {"input_data": {"dos": True}}}, "dos=False"),
        ({"matdyn_job": {"input_data": {"nk1": 10}}}, "DOS grid"),
        (
            {"matdyn_job": {"input_data": {"q_in_cryst_coord": False}}},
            "q_in_cryst_coord=True",
        ),
        ({"matdyn_job": {"input_data": {"q_in_band_form": True}}}, "explicit"),
        ({"matdyn_job": {"qpts": [[0, 0, 0]]}}, "flow's qpts"),
        ({"matdyn_job": {"input_data": {"l1": 2}}}, "SCF cell"),
        ({"matdyn_job": {"input_data": {"read_lr": True}}}, "together"),
        ({"q2r_job": {"input_data": {"write_lr": True}}}, "together"),
        ({"q2r_job": {"input_data": {"loto_2d": True}}}, "loto_2d"),
        ({"matdyn_job": {"input_data": {"loto_2d": True}}}, "loto_2d"),
    ],
)
def test_conflicting_parameters(tmp_path, params, message):
    with pytest.raises(ValueError, match=message):
        phonon_dispersion_flow(
            Atoms("H", cell=[2, 3, 4], pbc=True),
            prev_outdir=tmp_path,
            job_params=params,
        )


def test_requires_scf_save():
    with pytest.raises(ValueError, match="converged SCF"):
        phonon_dispersion_flow(Atoms("H", cell=[2, 3, 4], pbc=True))


@pytest.mark.parametrize("shape", [(2, 3), (1, 6)])
def test_mismatched_results(shape):
    path = {"matdyn_qpts": np.zeros((1, 3)), "natoms": 1}
    result = {
        "results": {
            "matdyn_results": {
                "frequencies": np.zeros(shape),
                "frequency_units": "cm^-1",
            }
        }
    }
    with pytest.raises(ValueError, match="dimensions"):
        summarize_phonon_dispersion(path, {}, {}, result)


def test_missing_frequencies():
    with pytest.raises(ValueError, match="did not produce"):
        summarize_phonon_dispersion({}, {}, {}, {"results": {}})


@pytest.mark.skipif(
    any(
        which(
            str(
                get_settings().ESPRESSO_BIN_DIR
                / get_settings().ESPRESSO_BINARIES[binary]
            )
        )
        is None
        for binary in ("pw", "ph", "q2r", "matdyn")
    ),
    reason="QE executables not configured",
)
def test_si_dispersion_qe(tmp_path, monkeypatch):
    """Exercise SCF -> 2x2x2 DFPT -> IFC -> automatic path with real QE.

    The settings test executable interoperability, not convergence.
    """
    from ase.build import bulk
    from monty.io import zopen

    monkeypatch.chdir(tmp_path)
    monkeypatch.setenv("OMP_NUM_THREADS", "1")
    copy_decompress_files(Path(__file__).parent / "data", "Si.upf.gz", tmp_path)
    atoms = bulk("Si", "diamond", a=5.43)
    with change_settings({"ESPRESSO_PSEUDO": tmp_path}):
        scf = static_job(
            atoms,
            preset=None,
            pseudopotentials={"Si": "Si.upf"},
            kpts=(4, 4, 4),
            input_data={
                "system": {"ecutwfc": 30, "ecutrho": 240, "occupations": "fixed"},
                "electrons": {"conv_thr": 1e-10},
            },
        )
        result = phonon_dispersion_flow(
            atoms,
            prev_outdir=scf["dir_name"],
            job_params={
                "phonon_job": {
                    "input_data": {"inputph": {"nq1": 2, "nq2": 2, "nq3": 2}}
                }
            },
        )
    data = result["dispersion"]
    path = prepare_phonon_path(atoms, True, 20, None, None, None)
    assert data["frequencies"].shape == (len(data["qpoints"]), 6)
    assert_array_equal(data["qpoints"], path["qpoints"])
    assert data["segments"] == path["segments"]
    assert data["metadata"]["qgrid"] == [2, 2, 2]
    # ph.x computes the dielectric tensor and Z* of insulating Si at Gamma.
    assert data["metadata"]["has_zstar"] is True
    assert data["metadata"]["loto_disable"] is False
    params = result["matdyn_job"]["parameters"]["input_data"]["input"]
    assert params["dos"] is False
    assert params["q_in_cryst_coord"] is True
    assert not params.get("q_in_band_form", False)
    assert not {"nk1", "nk2", "nk3"}.intersection(params)
    with zopen(
        Path(result["matdyn_job"]["dir_name"]) / "matdyn.freq.gp.gz", "rt"
    ) as handle:
        native = np.loadtxt(handle)
    assert len(native) == len(path["matdyn_qpts"])
    assert_allclose(
        data["frequencies"], native[path["matdyn_indices"], 1:], rtol=0, atol=1e-4
    )
    # Acoustic modes vanish at Gamma with the default crystal sum rule.
    gamma = [i for i, q in enumerate(data["qpoints"]) if not q.any()]
    assert gamma
    assert_allclose(data["frequencies"][gamma, :3], 0, atol=1)
    with zopen(Path(result["phonon_job"]["dir_name"]) / "ph.in.gz", "rt") as handle:
        _, cards = read_fortran_namelist(handle)
    assert cards == []

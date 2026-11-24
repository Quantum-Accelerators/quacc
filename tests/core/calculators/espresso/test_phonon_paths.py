from __future__ import annotations

from itertools import pairwise

import numpy as np
import pytest
from ase import Atoms
from ase.build import bulk
from ase.spacegroup import crystal
from numpy.testing import assert_allclose, assert_array_equal
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.symmetry.bandstructure import HighSymmKpath

from quacc.calculators.espresso.phonons import prepare_phonon_path


def matdyn_gamma_direction(rows, n):
    """Direction of the non-analytic term that matdyn.x uses at Gamma row n.

    This follows the q=0 branch of PHonon/PH/matdyn.f90 in QE 7.5 (0-based n).
    """
    if n == 0:
        direction = rows[0] - rows[1] if len(rows) > 1 else np.zeros(3)
    elif not rows[n - 1].any() and n < len(rows) - 1:
        direction = rows[n] - rows[n + 1]
    else:
        direction = rows[n] - rows[n - 1]
    norm = np.linalg.norm(direction)
    return direction / norm if norm else direction


def assert_gamma_directions(atoms, path):
    """Each Gamma point must use the direction of the branch it is plotted on."""
    reciprocal = 2 * np.pi * atoms.cell.reciprocal().array
    qpoints = path["qpoints"]
    rows = path["matdyn_qpts"] @ reciprocal
    assert_array_equal(path["matdyn_qpts"][path["matdyn_indices"]], qpoints)
    for start, stop in path["segments"]:
        for k in range(start, stop):
            if qpoints[k].any() or stop - start == 1:
                continue
            # Generated and well-formed manual paths have Gamma only at segment ends.
            neighbor = k + 1 if k == start else k - 1
            branch = (qpoints[neighbor] - qpoints[k]) @ reciprocal
            direction = matdyn_gamma_direction(rows, path["matdyn_indices"][k])
            assert abs(direction @ branch) == pytest.approx(np.linalg.norm(branch))


def test_continuous_endpoints():
    atoms = Atoms("H", cell=[2, 3, 4], pbc=True)
    points = [[0, 0, 0, 2], [0.5, 0, 0, 3], [0.5, 0.5, 0, 99]]
    result = prepare_phonon_path(atoms, False, 20, points, ["G", "X", "M"], None)
    expected = [
        [0, 0, 0],
        [0.25, 0, 0],
        [0.5, 0, 0],
        [0.5, 1 / 6, 0],
        [0.5, 1 / 3, 0],
        [0.5, 0.5, 0],
    ]
    assert_allclose(result["qpoints"], expected)
    # matdyn.x always receives the expanded, explicit list.
    assert_array_equal(result["matdyn_qpts"], result["qpoints"])
    assert result["matdyn_indices"] == list(range(6))
    assert_array_equal(result["qpoints"][[0, 2, 5]], np.asarray(points)[:, :3])
    assert result["segments"] == [(0, 6)]
    assert result["tick_indices"] == [0, 2, 5]
    assert result["labels"] == ["G", "", "X", "", "", "M"]
    assert_allclose(result["tick_positions"], [0, np.pi / 2, np.pi / 2 + np.pi / 3])


@pytest.mark.parametrize("endpoints", [False, True])
def test_disconnected_paths(endpoints):
    atoms = Atoms("H", cell=[2, 3, 4], pbc=True)
    if endpoints:
        points = [[0, 0, 0, 2], [0.5, 0, 0, 0], [0, 0.5, 0, 2], [0, 0.5, 0.5, 1]]
        labels = ["G", "X", "Y", "Z"]
        starts = None
    else:
        points = [
            [0, 0, 0],
            [0.25, 0, 0],
            [0.5, 0, 0],
            [0, 0.5, 0],
            [0, 0.5, 0.25],
            [0, 0.5, 0.5],
        ]
        labels = ["G", "", "X", "Y", "", "Z"]
        starts = [0, 3]
    result = prepare_phonon_path(atoms, False, 20, points, labels, starts)
    expected = [
        [0, 0, 0],
        [0.25, 0, 0],
        [0.5, 0, 0],
        [0, 0.5, 0],
        [0, 0.5, 0.25],
        [0, 0.5, 0.5],
    ]
    assert_allclose(result["qpoints"], expected)
    assert_allclose(result["matdyn_qpts"], expected)
    assert result["segments"] == [(0, 3), (3, 6)]
    assert result["tick_labels"] == ["G", "X", "Y", "Z"]
    assert_allclose(result["tick_positions"], [0, np.pi / 2, np.pi / 2, 3 * np.pi / 4])


def test_duplicate_gamma_for_directional_limits():
    atoms = Atoms("H", cell=[2, 3, 4], pbc=True)
    points = [[0.5, 0, 0, 2], [0, 0, 0, 1], [0, 0, 0, 2], [0, 0.5, 0, 1]]
    result = prepare_phonon_path(atoms, False, 20, points, ["X", "G", "G", "Y"], None)
    assert_array_equal(result["qpoints"][2:4], np.zeros((2, 3)))
    assert result["tick_indices"] == [0, 2, 3, 5]
    assert result["tick_positions"][1] == result["tick_positions"][2]
    # The repeated Gamma gives matdyn.x the incoming and then the outgoing limit.
    assert_array_equal(result["matdyn_qpts"], result["qpoints"])
    reciprocal = 2 * np.pi * atoms.cell.reciprocal().array
    rows = result["matdyn_qpts"] @ reciprocal
    assert_allclose(abs(matdyn_gamma_direction(rows, 2)), [1, 0, 0])
    assert_allclose(abs(matdyn_gamma_direction(rows, 3)), [0, 1, 0])


@pytest.mark.parametrize("endpoints", [False, True])
def test_gamma_starting_a_segment(endpoints):
    """After a break, Gamma must not take the LO-TO direction of the jump."""
    atoms = Atoms("H", cell=[2, 3, 4], pbc=True)
    if endpoints:
        points = [[0.5, 0, 0, 2], [0, 0.5, 0, 0], [0, 0, 0, 2], [0, 0, 0.5, 1]]
        labels, starts = ["X", "Y", "G", "Z"], None
    else:
        points = [[0.5, 0, 0], [0.25, 0.25, 0], [0, 0.5, 0], [0, 0, 0], [0, 0, 0.5]]
        labels, starts = ["X", "", "Y", "G", "Z"], [0, 3]
    result = prepare_phonon_path(atoms, False, 20, points, labels, starts)
    returned = result["qpoints"]
    assert result["segments"] == [(0, 3), (3, len(returned))]
    assert result["tick_labels"] == ["X", "Y", "G", "Z"]
    # One auxiliary Gamma precedes the segment and is omitted from the results.
    assert len(result["matdyn_qpts"]) == len(returned) + 1
    assert result["matdyn_indices"] == [0, 1, 2, *range(4, len(returned) + 1)]
    assert_array_equal(result["matdyn_qpts"][3], [0, 0, 0])
    assert_gamma_directions(atoms, result)


def test_gamma_ending_or_alone_in_a_segment():
    """Gamma at the end of a segment already uses the incoming direction."""
    atoms = Atoms("H", cell=[2, 3, 4], pbc=True)
    points = [[0.5, 0, 0], [0, 0, 0], [0, 0.5, 0], [0, 0, 0], [0, 0, 0.5]]
    result = prepare_phonon_path(atoms, False, 20, points, None, [0, 2, 3, 4])
    # Neither the segment-ending Gamma nor the single-point Gamma gets an extra row.
    assert_array_equal(result["matdyn_qpts"], points)
    assert result["matdyn_indices"] == [0, 1, 2, 3, 4]


def test_near_zero_coordinates_become_gamma():
    atoms = Atoms("H", cell=[2, 3, 4], pbc=True)
    points = [[0.5, 0, 0, 0], [1e-17, -2e-17, 0, 2], [0, 0.5, 0, 1]]
    result = prepare_phonon_path(atoms, False, 20, points, None, None)
    # matdyn.x tests the preceding q-point for exact zeros.
    assert not result["qpoints"][1].any()
    assert not result["matdyn_qpts"][1].any()
    assert_gamma_directions(atoms, result)


POLAR_UNIAXIAL = {
    "rutile": crystal(
        ["Ti", "O"],
        basis=[(0, 0, 0), (0.305, 0.305, 0)],
        spacegroup=136,
        cellpar=[4.59, 4.59, 2.96, 90, 90, 90],
    ),
    "hexagonal-BN": crystal(
        ["B", "N"],
        basis=[(1 / 3, 2 / 3, 0.25), (2 / 3, 1 / 3, 0.25)],
        spacegroup=194,
        cellpar=[2.5, 2.5, 6.66, 90, 90, 120],
    ),
    "wurtzite": bulk("GaN", "wurtzite", a=3.19, c=5.19),
}


@pytest.mark.parametrize("name", POLAR_UNIAXIAL)
def test_generated_gamma_lo_to_directions(name):
    """Anisotropic polar crystals: each Gamma limit follows its plotted branch."""
    atoms = POLAR_UNIAXIAL[name]
    result = prepare_phonon_path(atoms, True, 10, None, None, None)
    gamma = [i for i, point in enumerate(result["qpoints"]) if not point.any()]
    assert gamma
    for i in gamma:
        assert result["labels"][i] == "Γ"
    assert_gamma_directions(atoms, result)


@pytest.mark.parametrize(
    "atoms",
    [
        bulk("Si", "diamond", a=5.43),
        bulk("Al", "fcc", a=4, cubic=True),
        bulk("Mg", "hcp", a=3.2, c=5.2),
        Atoms("H", cell=[5, 4, 3.2], pbc=True),
        Atoms("H", cell=[[3, 0, 0], [0.7, 4, 0], [0.4, 0.2, 5]], pbc=True),
        Atoms("H", cell=[[3, 0, 0], [6, 4, 0], [0, 0, 5]], pbc=True),
    ],
    ids=[
        "diamond",
        "conventional-fcc",
        "hexagonal",
        "orthorhombic",
        "triclinic",
        "unreduced",
    ],
)
def test_auto_path_original_basis(atoms):
    atoms = atoms.copy()
    atoms.rotate(27, "z", rotate_cell=True)
    original = atoms.copy()
    line_density = 2
    result = prepare_phonon_path(atoms, True, line_density, None, None, None)
    reference_structure = AseAtomsAdaptor.get_structure(atoms).get_reduced_structure(
        "LLL"
    )
    reference = HighSymmKpath(reference_structure, path_type="latimer_munro")
    reciprocal = 2 * np.pi * atoms.cell.reciprocal().array
    qpoints = result["qpoints"]

    def cartesian(label):
        return reference.rec_lattice.get_cartesian_coords(
            reference.kpath["kpoints"][label]
        )

    for index, label in zip(result["tick_indices"], result["tick_labels"], strict=True):
        assert_allclose(qpoints[index] @ reciprocal, cartesian(label), atol=1e-10)

    # Interior Gamma points are repeated at the same distance; otherwise the
    # labels and point counts follow pymatgen's branches and line density.
    ticks = list(zip(result["tick_labels"], result["tick_positions"], strict=True))
    collapsed = [
        label
        for i, (label, position) in enumerate(ticks)
        if not (i and label == ticks[i - 1][0] and position == ticks[i - 1][1])
    ]
    expected_labels = [
        str(label) for branch in reference.kpath["path"] for label in branch
    ]
    assert collapsed == expected_labels
    interior_gamma = sum(
        str(label) == "Γ"
        for branch in reference.kpath["path"]
        for label in branch[1:-1]
    )
    npoints = interior_gamma
    for branch in reference.kpath["path"]:
        npoints += 1
        for start, end in pairwise(branch):
            distance = np.linalg.norm(cartesian(end) - cartesian(start))
            npoints += max(1, int(np.ceil(distance * line_density)))
    assert len(qpoints) == npoints
    assert len(result["segments"]) == len(reference.kpath["path"]) + interior_gamma
    for start, stop in result["segments"]:
        increments = np.linalg.norm(
            np.diff(qpoints[start:stop], axis=0) @ reciprocal, axis=1
        )
        assert_allclose(np.diff(result["distances"][start:stop]), increments)
        if start:
            assert result["distances"][start] == result["distances"][start - 1]
    assert_gamma_directions(atoms, result)
    assert result["natoms"] == len(atoms)
    assert all(isinstance(index, int) for pair in result["segments"] for index in pair)
    assert_allclose(result["cell"], original.cell)
    assert_allclose(atoms.positions, original.positions)
    assert_allclose(atoms.cell, original.cell)


def test_manual_path_for_2d_cell():
    atoms = Atoms(
        "BN",
        scaled_positions=[[0, 0, 0.5], [1 / 3, 2 / 3, 0.5]],
        cell=[[2.5, 0, 0], [-1.25, 2.165, 0], [0, 0, 15]],
        pbc=[True, True, False],
    )
    points = [[0, 0, 0, 3], [0.5, 0, 0, 2], [1 / 3, 1 / 3, 0, 4], [0, 0, 0, 1]]
    result = prepare_phonon_path(atoms, False, 20, points, ["G", "M", "K", "G"], None)
    assert not result["qpoints"][:, 2].any()
    with pytest.raises(ValueError, match="in-plane path"):
        prepare_phonon_path(atoms, True, 20, None, None, None)


def test_numpy_segment_starts():
    result = prepare_phonon_path(
        Atoms("H", cell=[2, 3, 4], pbc=True),
        False,
        20,
        np.array([[0, 0, 0], [0.5, 0, 0], [0, 0.5, 0]]),
        None,
        np.array([0, 2]),
    )
    assert result["segments"] == [(0, 2), (2, 3)]
    assert all(type(index) is int for pair in result["segments"] for index in pair)


@pytest.mark.parametrize(
    ("points", "labels", "starts", "message"),
    [
        (None, None, None, "Supply qpts"),
        ([], None, None, "shape"),
        (0, None, None, "shape"),
        ([[0, 0, 0]], ["G", "X"], None, "one string"),
        ([[0, 0, 0]], [0], None, "one string"),
        ([[0, 0, 0]], None, [1], "beginning with 0"),
        ([[0, 0, 0]], None, [], "beginning with 0"),
        ([[0, 0, 0]], None, [0, 0], "increasing"),
        ([[0, 0, 0]], None, [0, 1], "indices"),
        ([[0, 0, 0]], None, [0.0], "indices"),
        ([[0, 0, 0]], None, [False], "indices"),
        ([[0, 0, 0, 2]], None, [0], "zero interpolation counts"),
        ([[0, 0, 0, -1], [0.5, 0, 0, 1]], None, None, "nonnegative integers"),
        ([[0, np.nan, 0]], None, None, "finite"),
    ],
)
def test_invalid_manual_path(points, labels, starts, message):
    with pytest.raises(ValueError, match=message):
        prepare_phonon_path(
            Atoms("H", cell=[2, 3, 4], pbc=True), False, 20, points, labels, starts
        )


@pytest.mark.parametrize("density", [0, -1, np.nan])
def test_invalid_density(density):
    with pytest.raises(ValueError, match="line_density"):
        prepare_phonon_path(bulk("Al"), True, density, None, None, None)

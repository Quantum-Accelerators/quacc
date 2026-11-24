"""Path preparation and result assembly for interpolated DFPT dispersions."""

from __future__ import annotations

from itertools import pairwise
from numbers import Integral
from typing import TYPE_CHECKING, cast

import numpy as np
from ase.io.espresso import Namelist
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.symmetry.bandstructure import HighSymmKpath

from quacc.calculators.espresso.utils import matdyn_qpoint_cards
from quacc.utils.dicts import recursive_dict_merge

if TYPE_CHECKING:
    from collections.abc import Sequence
    from typing import Any

    from ase.atoms import Atoms
    from numpy.typing import ArrayLike, NDArray
    from pymatgen.core.structure import Structure

    from quacc.types import EspressoPhononDispersionSchema, RunSchema


def prepare_phonon_path(
    atoms: Atoms,
    make_bandpath: bool,
    line_density: float,
    qpts: ArrayLike | None,
    qpoint_labels: Sequence[str] | None,
    segment_starts: Sequence[int] | None,
) -> dict[str, Any]:
    """Prepare fractional q-points, the matdyn.x q-point list, and plot metadata.

    Endpoint paths are expanded with QE's interval counts, and matdyn.x always
    receives an explicit list. Segment ranges use Python's [start, stop)
    convention, and distances exclude jumps between segments.

    At Gamma, matdyn.x takes the direction of the non-analytic (LO-TO) term from
    the preceding q-point unless that point is also Gamma, in which case it uses
    the following point. Where Gamma starts a segment after a break, an auxiliary
    Gamma point is placed before it in the matdyn.x list, so the outgoing direction
    is used; its frequencies are discarded through ``matdyn_indices``. Generated
    paths also split at interior Gamma points to keep both limiting directions.
    """
    if len(atoms) == 0 or atoms.cell.rank != 3:
        raise ValueError(
            "Phonon dispersion requires a nonempty cell with three lattice vectors."
        )
    reciprocal = 2 * np.pi * atoms.cell.reciprocal().array
    if qpts is None:
        if not make_bandpath:
            raise ValueError("Supply qpts when make_bandpath=False.")
        if qpoint_labels is not None or segment_starts is not None:
            raise ValueError("Manual labels and segment_starts require qpts.")
        points, labels = _latimer_munro_endpoints(atoms, line_density, reciprocal)
        source = "latimer_munro"
    else:
        points = np.asarray(qpts, dtype=float)
        if points.ndim != 2 or points.shape[1] not in (3, 4):
            raise ValueError("qpts must have shape (n, 3) or (n, 4).")
        labels = (
            list(qpoint_labels) if qpoint_labels is not None else [""] * len(points)
        )
        source = "manual"
    band_form = points.shape[1] == 4
    matdyn_qpoint_cards(points, q_in_band_form=band_form)
    if len(labels) != len(points) or any(
        not isinstance(label, str) for label in labels
    ):
        raise ValueError(
            "qpoint_labels must contain one string for each input q-point."
        )

    if band_form:
        if segment_starts is not None:
            raise ValueError(
                "For endpoint paths, use zero interpolation counts for path breaks."
            )
        qpoints, labels, starts = _expand_endpoints(points, labels)
    else:
        qpoints = points.copy()
        starts = _check_segment_starts(segment_starts, len(points))
    # matdyn.x recognizes a preceding Gamma point only if it is exactly zero.
    qpoints[np.abs(qpoints) < 1e-10] = 0.0
    # TODO: Manual paths are not split at Gamma. An opt-in split would give both
    # LO limits without requiring users to repeat Gamma by hand.
    if source == "latimer_munro":
        qpoints, labels, starts = _split_at_gamma(qpoints, labels, starts)
    matdyn_qpts, matdyn_indices = _matdyn_qpoints(qpoints, starts)

    segments = list(zip(starts, [*starts[1:], len(qpoints)], strict=True))
    increments = np.linalg.norm(np.diff(qpoints, axis=0) @ reciprocal, axis=1)
    increments[np.asarray(starts[1:], dtype=int) - 1] = 0
    distances = np.r_[0.0, np.cumsum(increments)]
    tick_indices = [i for i, label in enumerate(labels) if label]
    return {
        "qpoints": qpoints,
        "matdyn_qpts": matdyn_qpts,
        "matdyn_indices": matdyn_indices,
        "distances": distances,
        "labels": labels,
        "tick_indices": tick_indices,
        "tick_labels": [labels[i] for i in tick_indices],
        "tick_positions": distances[tick_indices],
        "segments": segments,
        "cell": atoms.cell.array.copy(),
        "natoms": len(atoms),
        "path_source": source,
    }


def _latimer_munro_endpoints(
    atoms: Atoms, line_density: float, reciprocal: NDArray
) -> tuple[NDArray, list[str]]:
    """Generate Latimer-Munro endpoints and counts in the reciprocal basis of atoms."""
    # TODO: Only Latimer-Munro is offered, as in bands_flow. Other conventions
    # (e.g. setyawan_curtarolo, hinuma) and automatic in-plane paths for 2D
    # cells would need a path_type option shared with quacc.utils.kpts.
    if not np.all(atoms.pbc):
        raise ValueError(
            "Automatic paths require a fully periodic cell. For 2D systems, "
            "supply an in-plane path through qpts."
        )
    if not np.isfinite(line_density) or line_density <= 0:
        raise ValueError("line_density must be finite and positive.")
    structure = cast("Structure", AseAtomsAdaptor.get_structure(atoms))
    # LLL only changes the basis of this cell. It avoids Latimer-Munro's
    # warning about unreduced lattices without converting to a primitive cell.
    reduced = structure.get_reduced_structure(reduction_algo="LLL")
    transform = reduced.lattice.matrix @ np.linalg.inv(atoms.cell.array)
    if (
        len(reduced) != len(atoms)
        or not np.allclose(transform, np.rint(transform), rtol=0, atol=1e-8)
        or not np.isclose(abs(np.linalg.det(transform)), 1, rtol=0, atol=1e-8)
    ):
        raise ValueError(
            "Path preparation must preserve the SCF cell volume and atom count."
        )
    kpath = HighSymmKpath(
        reduced,
        path_type="latimer_munro",
        has_magmoms=bool(np.any(structure.site_properties.get("magmom", 0))),
    )
    if not np.allclose(
        kpath.rec_lattice.matrix,
        reduced.lattice.reciprocal_lattice.matrix,
        rtol=1e-8,
        atol=1e-8,
    ):
        raise ValueError("The generated path has an inconsistent reciprocal basis.")
    special_points = {}
    for label, point in kpath.kpath["kpoints"].items():
        cartesian = kpath.rec_lattice.get_cartesian_coords(point)
        fractional = cartesian @ np.linalg.inv(reciprocal)
        if not np.allclose(fractional @ reciprocal, cartesian, rtol=1e-8, atol=1e-8):
            raise ValueError(
                "Could not transform the q-path to the SCF reciprocal basis."
            )
        special_points[label] = fractional
    endpoints = []
    labels = []
    for branch in kpath.kpath["path"]:
        for i, label in enumerate(branch):
            point = special_points[label]
            count = 0
            if i + 1 < len(branch):
                following = special_points[branch[i + 1]]
                distance = np.linalg.norm((following - point) @ reciprocal)
                count = max(1, int(np.ceil(distance * line_density)))
            endpoints.append([*point, count])
            labels.append(str(label))
    return np.asarray(endpoints), labels


def _expand_endpoints(
    points: NDArray, labels: list[str]
) -> tuple[NDArray, list[str], list[int]]:
    """Expand endpoints like QE: n intervals to the next point, or a break at 0."""
    qpoints = [points[0, :3]]
    expanded_labels = [labels[0]]
    starts = [0]
    for i, (start, end) in enumerate(pairwise(points)):
        count = int(start[3])
        if count == 0:
            starts.append(len(qpoints))
        else:
            qpoints.extend(
                start[:3] + j / count * (end[:3] - start[:3]) for j in range(1, count)
            )
            expanded_labels.extend([""] * (count - 1))
        qpoints.append(end[:3])
        expanded_labels.append(labels[i + 1])
    return np.asarray(qpoints), expanded_labels, starts


def _check_segment_starts(
    segment_starts: Sequence[int] | None, npoints: int
) -> list[int]:
    """Validate the first index of each explicit segment."""
    starts = list(segment_starts) if segment_starts is not None else [0]
    if (
        not starts
        or any(
            not isinstance(index, Integral) or isinstance(index, bool)
            for index in starts
        )
        or starts[0] != 0
        or any(index < 0 or index >= npoints for index in starts)
        or any(b <= a for a, b in pairwise(starts))
    ):
        raise ValueError(
            "segment_starts must be increasing point indices beginning with 0."
        )
    return [int(index) for index in starts]


def _segment_stops(starts: list[int], npoints: int) -> set[int]:
    """Return the index of the last point in each segment."""
    return {start - 1 for start in starts[1:]} | {npoints - 1}


def _split_at_gamma(
    qpoints: NDArray, labels: list[str], starts: list[int]
) -> tuple[NDArray, list[str], list[int]]:
    """Repeat interior Gamma points so that each side of Gamma has its own limit."""
    breaks = set(starts)
    stops = _segment_stops(starts, len(qpoints))
    rows: list[NDArray] = []
    new_labels: list[str] = []
    new_starts: list[int] = []
    for i, (point, label) in enumerate(zip(qpoints, labels, strict=True)):
        if i in breaks:
            new_starts.append(len(rows))
        elif i not in stops and not point.any() and qpoints[i - 1].any():
            # This copy ends the incoming segment; the next starts the outgoing one.
            rows.append(point)
            new_labels.append(label)
            new_starts.append(len(rows))
        rows.append(point)
        new_labels.append(label)
    return np.asarray(rows), new_labels, new_starts


def _matdyn_qpoints(qpoints: NDArray, starts: list[int]) -> tuple[NDArray, list[int]]:
    """Insert an auxiliary Gamma before each Gamma that starts a continuing segment."""
    # This follows the Gamma direction rule of matdyn.f90 in QE 7.5. If a later
    # QE release changes that rule, revisit this function and _split_at_gamma.
    breaks = set(starts[1:])
    stops = _segment_stops(starts, len(qpoints))
    rows: list[NDArray] = []
    indices = []
    for i, point in enumerate(qpoints):
        if i in breaks and i not in stops and not point.any() and qpoints[i - 1].any():
            rows.append(np.zeros(3))
        indices.append(len(rows))
        rows.append(point)
    return np.asarray(rows), indices


def validate_phonon_dispersion_params(
    defaults: dict[str, dict[str, Any]], job_params: dict[str, Any] | None
) -> None:
    """Validate the effective ph.x, q2r.x, and matdyn.x inputs of the flow."""
    params = {}
    swaps = job_params or {}
    for name, binary in [
        ("phonon_job", "ph"),
        ("q2r_job", "q2r"),
        ("matdyn_job", "matdyn"),
    ]:
        params[name] = recursive_dict_merge(
            defaults.get(name), swaps.get("all"), swaps.get(name)
        )
        inputs = Namelist(params[name].get("input_data"))
        inputs.to_nested(binary=binary)
        params[name]["input_data"] = inputs
    ph = params["phonon_job"]["input_data"]["inputph"]
    matdyn = params["matdyn_job"]["input_data"]["input"]
    q2r = params["q2r_job"]["input_data"]["input"]
    if (
        not ph.get("ldisp", False)
        or ph.get("qplot", False)
        or params["phonon_job"].get("use_phcg", False)
    ):
        raise ValueError(
            "Phonon dispersion requires ph.x with ldisp=True and qplot=False."
        )
    if any(
        not isinstance(ph.get(key), Integral)
        or isinstance(ph[key], bool)
        or ph[key] <= 0
        for key in ("nq1", "nq2", "nq3")
    ):
        raise ValueError(
            "The coarse DFPT nq1/nq2/nq3 grid must contain positive integers."
        )
    if matdyn.get("dos", False) or not matdyn.get("q_in_cryst_coord", False):
        raise ValueError(
            "Phonon dispersion requires dos=False and q_in_cryst_coord=True."
        )
    if any(key in matdyn for key in ("nk1", "nk2", "nk3")):
        raise ValueError("nk1/nk2/nk3 describe the DOS grid, not a dispersion path.")
    if any(key in params["matdyn_job"] for key in ("qpts", "additional_cards")) or (
        matdyn.get("q_in_band_form", False)
    ):
        raise ValueError(
            "Set the dispersion path through the flow's qpts argument; "
            "matdyn.x receives an explicit q-point list."
        )
    if (
        matdyn.get("readtau", False)
        or any(matdyn.get(key, 1) != 1 for key in ("l1", "l2", "l3"))
        or any(key.startswith("at(") or key == "at" for key in matdyn)
    ):
        raise ValueError(
            "The matdyn cell must match the SCF cell; supercell overrides are unsupported."
        )
    # A text force-constant file has an extra long-range column when write_lr=True,
    # and matdyn.x reads that column only with read_lr=True.
    if bool(q2r.get("write_lr", False)) != bool(matdyn.get("read_lr", False)):
        raise ValueError(
            "Set q2r write_lr=True and matdyn read_lr=True together for long-range IFCs."
        )
    # q2r.x removes, and matdyn.x restores, the dipole-dipole term with loto_2d.
    if bool(q2r.get("loto_2d", False)) != bool(matdyn.get("loto_2d", False)):
        raise ValueError("Set loto_2d identically for q2r_job and matdyn_job.")


def summarize_phonon_dispersion(
    path: dict[str, Any],
    phonon_results: RunSchema,
    q2r_results: RunSchema,
    matdyn_results: RunSchema,
) -> EspressoPhononDispersionSchema:
    """Assemble resolved results without evaluating workflow futures in the flow."""
    parsed = matdyn_results["results"].get("matdyn_results", {})
    if "frequencies" not in parsed:
        raise ValueError("matdyn did not produce parsed phonon frequencies.")
    frequencies = np.asarray(parsed["frequencies"])
    if frequencies.shape != (len(path["matdyn_qpts"]), 3 * path["natoms"]):
        raise ValueError(
            "matdyn frequency dimensions do not match the requested path and SCF atom count."
        )
    if parsed.get("frequency_units") != "cm^-1":
        raise ValueError("Expected matdyn frequencies in cm^-1.")
    if "qpoints" in parsed:
        _check_matdyn_qpoints(path, np.asarray(parsed["qpoints"], dtype=float))
    ph_input = phonon_results["parameters"]["input_data"]["inputph"]
    q2r_input = q2r_results["parameters"]["input_data"]["input"]
    matdyn_input = matdyn_results["parameters"]["input_data"]["input"]
    # TODO: Frequencies keep QE's ascending order at each q-point; connecting
    # branches across crossings would require parsing the flvec eigenvectors.
    return {
        "phonon_job": phonon_results,
        "q2r_job": q2r_results,
        "matdyn_job": matdyn_results,
        "dispersion": {
            "qpoints": path["qpoints"],
            "qpoint_coordinate_system": "fractional",
            "frequencies": frequencies[path["matdyn_indices"]],
            "frequency_units": "cm^-1",
            "distances": path["distances"],
            "distance_units": "angstrom^-1",
            "labels": path["labels"],
            "tick_indices": path["tick_indices"],
            "tick_labels": path["tick_labels"],
            "tick_positions": path["tick_positions"],
            "segments": path["segments"],
            "metadata": {
                "cell": path["cell"],
                "natoms": path["natoms"],
                "path_source": path["path_source"],
                "qgrid": [ph_input[key] for key in ("nq1", "nq2", "nq3")],
                "zasr": q2r_input.get("zasr", "no"),
                "asr": matdyn_input.get("asr", "no"),
                "write_lr": q2r_input.get("write_lr", False),
                "read_lr": matdyn_input.get("read_lr", False),
                "loto_disable": matdyn_input.get("loto_disable", False),
                "loto_2d": matdyn_input.get("loto_2d", False),
                "has_zstar": parsed.get("has_zstar"),
            },
        },
    }


def _check_matdyn_qpoints(path: dict[str, Any], cartesian: NDArray) -> None:
    """Check matdyn's Cartesian q-points against the path in the cell of ``atoms``.

    matdyn.x writes q-points in units of 2pi/alat and converts fractional input
    with the cell of the force-constant file. Agreement up to one positive scale
    factor confirms that the path and the SCF calculation share a reciprocal basis.
    """
    expected = path["matdyn_qpts"] @ (2 * np.pi * np.linalg.inv(path["cell"]).T)
    if cartesian.shape != expected.shape:
        raise ValueError("matdyn q-points do not match the requested path.")
    norm = np.vdot(cartesian, cartesian)
    scale = np.vdot(cartesian, expected) / norm if norm else 0.0
    if (norm and scale <= 0) or not np.allclose(
        scale * cartesian, expected, rtol=0, atol=1e-4
    ):
        raise ValueError(
            "matdyn q-points are inconsistent with the cell of atoms. "
            "Use the same cell and orientation as the SCF calculation."
        )

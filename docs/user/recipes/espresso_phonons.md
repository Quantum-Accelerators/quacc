# Quantum ESPRESSO phonon dispersions

[quacc.recipes.espresso.phonons.phonon_dispersion_flow][] interpolates harmonic DFPT phonons using `ph.x`, `q2r.x`, and `matdyn.x`. It reads a converged `pw.x` SCF save and retains the SCF cell and atom ordering. Configure the QE binaries and pseudopotential directory through the usual [Quacc settings](../settings/settings.md).

## Silicon example

The following example uses a two-atom diamond primitive cell and a PBEsol norm-conserving `Si.upf`. Phonons are computed at the equilibrium structure, so relax the structure with the same pseudopotential and cutoffs before the SCF. The cutoffs and grids are an execution example; converge them for the chosen pseudopotential and material.

```python
from ase.build import bulk

from quacc.recipes.espresso.core import static_job
from quacc.recipes.espresso.phonons import phonon_dispersion_flow

atoms = bulk("Si", "diamond", a=5.4437023729394526)
scf = static_job(
    atoms,
    preset=None,
    pseudopotentials={"Si": "Si.upf"},
    kpts=(8, 8, 8),
    input_data={
        "system": {"ecutwfc": 40, "occupations": "fixed"},
        "electrons": {"conv_thr": 1e-12},
    },
)
result = phonon_dispersion_flow(
    atoms,
    prev_outdir=scf["dir_name"],
    line_density=20,
    job_params={
        "phonon_job": {"input_data": {"inputph": {"nq1": 3, "nq2": 3, "nq3": 3}}}
    },
)
dispersion = result["dispersion"]
```

`occupations="fixed"` treats Si as an insulator; Quacc then logs that the default smearing settings of `static_job` are ignored. Without `job_params`, the DFPT grid is 4×4×4. The 3×3×3 grid above is not converged: X and L are not on it, and with the 4×4×4 grid the transverse acoustic frequency at X rises from 119 to 137 cm⁻¹. Compare dispersions from two grids before relying on one.

The four steps have distinct sampling inputs:

| Program | Sampling | Role |
| --- | --- | --- |
| `pw.x` | Electronic `kpts` | Converged SCF save |
| `ph.x` | Uniform `nq1/nq2/nq3`, `ldisp=True` | Coarse-grid dynamical matrices |
| `q2r.x` | The complete dynamical-matrix set | Real-space force constants |
| `matdyn.x` | High-symmetry path, `dos=False` | Interpolated frequencies |

Pass the same `atoms` as in the SCF, with the same cell and orientation. After `matdyn.x` runs, the flow compares the Cartesian q-points it wrote with the requested path and raises an error if the reciprocal bases do not match. Instead of `prev_outdir`, the SCF files can be supplied through `copy_files`. Customize `phonon_job`, `q2r_job`, and `matdyn_job` through `job_params` and `job_decorators`, as in the other Espresso flows. Path preparation and result assembly also run as jobs, so that a workflow engine can pass deferred atoms and results to them; `job_decorators` does not apply to these two jobs. The flow defines the path input of `matdyn.x` itself (`dos=False`, `q_in_cryst_coord=True`, and an explicit q-point list). Before any calculation runs, it rejects settings that conflict with a dispersion, such as `ldisp=False`, a DOS grid (`nk1/nk2/nk3`), a path set through `matdyn_job` instead of `qpts`, or a `matdyn.x` supercell.

## Automatic paths

Without `qpts`, the flow generates pymatgen's Latimer–Munro path, as `bands_flow()` does, and returns its labels. These are not the textbook letters: for the Si example, `d`, `a` and `r` are X, L and W, and `q_{1}` and `q` are the equivalent points K and U. For textbook labels, supply a [manual path](#manual-paths). The path is generated on a reduced basis of the same lattice (same volume and atom count), and its q-points are transformed back into the reciprocal basis of the SCF cell; the structure is not converted to a standard or primitive cell. Automatic paths require a fully periodic cell; see [Two-dimensional materials](#two-dimensional-materials) for layers in a cell with vacuum.

`line_density` sets the number of interpolation intervals per inverse angstrom (including `2π`) on the automatic path. A higher value gives more points along the path, but the frequencies are only as accurate as the force constants from the DFPT grid.

## Manual paths

Supply fractional reciprocal coordinates in the SCF basis through `qpts`, in one of the two forms below. Supplying `qpts` takes precedence over automatic path generation; `make_bandpath=False` requires a manual path. In both forms, `qpoint_labels` optionally gives one string per row of `qpts`; use empty strings for unlabeled rows. Labeled rows become the axis ticks.

**Endpoints with interpolation counts (four columns).** Each row is an endpoint followed by a QE interpolation count. A count of `n` adds `n` intervals to the following endpoint, and the last count is ignored. For example, a continuous five-endpoint path with counts `[20, 20, 20, 20, 1]` has 81 output points. A zero count marks a break: both endpoints are kept, and nothing is interpolated across the gap. In this form, breaks are set only by zero counts; `segment_starts` is rejected. The flow expands the endpoints before `matdyn.x` runs, and the interpolated points have empty labels.

ASE gives the special points of a cell in its own reciprocal basis. Take all points from one call: symmetry-equivalent copies chosen by hand need not be neighbors on the zone boundary. A textbook path for the Si example:

```python
special = atoms.cell.bandpath("GXWKGLUWLK,UX", npoints=0).special_points
route = ["G", "X", "W", "K", "G", "L", "U", "W", "L", "K", "U", "X"]
counts = [20, 10, 7, 21, 17, 13, 7, 14, 13, 0, 7, 1]  # 0: break between K and U
endpoints = [[*special[point], count] for point, count in zip(route, counts)]
labels = ["Γ" if point == "G" else point for point in route]
result = phonon_dispersion_flow(
    atoms, prev_outdir=scf["dir_name"], qpts=endpoints, qpoint_labels=labels
)
```

**Explicit points (three columns).** Each row is evaluated as given, and the flow adds no points. Without `segment_starts`, all points form one connected segment, so consecutive points are joined in plots and counted in `distances`. To break the path, pass `segment_starts`: the row index at which each segment begins, in increasing order and starting with 0. For example, Γ–X and L–Γ as two segments with 20 intervals each:

```python
import numpy as np

gamma_x = np.linspace(special["G"], special["X"], 21)
l_gamma = np.linspace(special["L"], special["G"], 21)
points = np.vstack([gamma_x, l_gamma])
labels = ["Γ", *[""] * 19, "X", "L", *[""] * 19, "Γ"]
result = phonon_dispersion_flow(
    atoms,
    prev_outdir=scf["dir_name"],
    qpts=points,
    qpoint_labels=labels,
    segment_starts=[0, 21],
)
```

**Gamma in manual paths.** In polar insulators, the LO frequencies at Gamma depend on the direction from which q approaches Gamma (see [Sum rules and polar materials](#sum-rules-and-polar-materials)). The flow does not split manual paths at Gamma, so a single interior Gamma point gives only the limit from the incoming direction. To get both limits, repeat Gamma and start a new segment at the second copy. In the endpoint form, give the first copy a count of 0, for example `[..., [0, 0, 0, 0], [0, 0, 0, 20], ...]`; a nonzero count would instead interpolate between the two copies and add more Gamma points. In the explicit form, list Gamma twice and add the index of the second copy to `segment_starts`. Materials without LO-TO splitting, such as Si, do not need this.

To evaluate another path without repeating `ph.x`, pass the force constants of a finished run to `matdyn_job`; the example below reuses the endpoint path above. It returns QE's raw output (Cartesian q-points in `2π/alat`, without labels or segments), does not apply the Gamma handling described below, and needs `asr` set explicitly:

```python
from quacc.recipes.espresso.phonons import matdyn_job

matdyn = matdyn_job(
    [{"source": result["q2r_job"]["dir_name"], "filenames": "q2r.fc*"}],
    qpts=endpoints,
    input_data={
        "input": {
            "asr": "crystal",
            "dos": False,
            "q_in_band_form": True,
            "q_in_cryst_coord": True,
        }
    },
)
frequencies = matdyn["results"]["matdyn_results"]["frequencies"]
```

## Structured results

The return value follows [quacc.types.EspressoPhononDispersionSchema][] and includes `phonon_job`, `q2r_job`, `matdyn_job`, and `dispersion`.

| Field in `dispersion` | Meaning |
| --- | --- |
| `qpoints` | Fractional reciprocal coordinates, shape `(nqpoints, 3)` |
| `frequencies` | Frequencies in cm⁻¹, shape `(nqpoints, 3 * natoms)`; negative modes retained |
| `distances` | Cumulative reciprocal-space distance in Å⁻¹, including `2π` and omitting inter-segment jumps |
| `labels` | One string per output point; empty strings between labeled endpoints |
| `tick_indices`, `tick_positions`, `tick_labels` | Endpoint indices, distances, and labels for axis ticks |
| `segments` | Half-open `[start, stop)` ranges for connected path components; generated paths also split at interior Gamma points |
| `metadata` | SCF cell, atom count, path source, coarse grid, ASR, long-range IFC and LO-TO settings, and `has_zstar` |

Use each segment's range separately when consuming the data, and use `tick_positions` with `tick_labels`. At a break, two distinct endpoint labels can share the same distance, and a repeated interior Gamma point gives two identical ticks. The following matplotlib example joins the labels at a shared position (`d|Γ`) and wraps them in `$...$`, so that Latimer–Munro subscripts such as `q_{1}` are typeset:

```python
import matplotlib.pyplot as plt

fig, ax = plt.subplots()
x, freqs = dispersion["distances"], dispersion["frequencies"]
for start, stop in dispersion["segments"]:
    ax.plot(x[start:stop], freqs[start:stop], color="tab:blue")
ticks = {}
for position, label in zip(dispersion["tick_positions"], dispersion["tick_labels"]):
    ticks.setdefault(position, []).append(label)
names = ["|".join(dict.fromkeys(labels)) for labels in ticks.values()]
ax.set_xticks(list(ticks), [f"${name}$" for name in names])
ax.set_ylabel("Frequency (cm$^{-1}$)")
fig.savefig("phonons.png")
```

The existing `matdyn_job` result also retains QE's raw output for every q-point it evaluated, including any auxiliary Gamma point described below, with Cartesian q-points and `.freq.gp` distances in `2π/alat` units. QE's `.gp` distances include inter-segment jumps, so they differ from the flow's distances. Frequencies are read from `.freq`; the companion `.freq.gp`, when present, is cross-checked without reordering modes. See the official [matdyn input specification](https://www.quantum-espresso.org/Doc/INPUT_MATDYN.html).

For JSON export, apply `monty.json.jsanitize` to the combined result before encoding it. The preserved `ph.x` job dictionary can contain tuple keys.

## Sum rules and polar materials

The flow uses `zasr="crystal"` in `q2r.x` and `asr="crystal"` in `matdyn.x`, matching `phonon_dos_flow()`. Override these through `job_params`. No dimensional correction is enabled automatically.

In a dispersion calculation, `ph.x` computes the dielectric tensor and Born effective charges at Gamma automatically when no smearing or tetrahedra are used, so `matdyn.x` applies the non-analytic (LO-TO) term for insulators. Choose `loto_disable` or `loto_2d` explicitly where appropriate; `loto_2d` must be set identically for `q2r_job` and `matdyn_job`, because `q2r.x` removes and `matdyn.x` restores the dipole-dipole term with it. For long-range IFCs, supply `write_lr=True` to `q2r_job` and `read_lr=True` to `matdyn_job` together, since `write_lr` adds a column to the force-constant file. These are required for `asr="all"` in infrared-active solids; see the [q2r](https://www.quantum-espresso.org/Doc/INPUT_Q2R.html) and [matdyn](https://www.quantum-espresso.org/Doc/INPUT_MATDYN.html) specifications.

`metadata` records `loto_disable` and `loto_2d` as passed to `matdyn.x`, and `has_zstar`, which is true when the force constants contain the dielectric tensor and Born effective charges. `has_zstar` is read from text force-constant files and is `None` otherwise. `matdyn.x` applies the LO-TO term at Gamma only when `has_zstar` is true and neither `loto_disable` nor `loto_2d` is set. `loto_disable` changes only the Gamma rows; elsewhere, the dipole-dipole term is still included whenever `has_zstar` is true.

At Gamma, `matdyn.x` takes the direction of the non-analytic term from the preceding q-point, or from the following one when the preceding point is also Gamma. In anisotropic polar crystals the LO frequencies at Gamma depend on this direction. Generated paths therefore repeat every interior Gamma point: one copy ends the incoming segment and the next starts the outgoing segment, at the same distance. Where Gamma starts a segment after a break, the flow evaluates an auxiliary Gamma point first and omits it from `dispersion`, so the outgoing direction is used. Manual paths are evaluated as given, except for this auxiliary point; to obtain both limits, repeat Gamma as described in [Manual paths](#manual-paths).

## Two-dimensional materials

QE treats a layer in a cell with vacuum as periodic in all three directions. Set `atoms.pbc = True` to use the automatic path. The path then also contains legs with an out-of-plane component `q_z`. With `nq3=1`, the force constants contain no interlayer terms, so these legs repeat in-plane frequencies exactly; the vacuum must still be thick enough that periodic images of the layer do not interact in the SCF and DFPT calculations. The legs with `q_z = 0` form the in-plane loop Gamma–M–K–Gamma (conventional labels); alternatively, supply an in-plane path through `qpts`.

The out-of-plane acoustic (flexural, ZA) branch of a 2D material is quadratic near Gamma, `ω_ZA ∝ |q|²`. This follows from rotational invariance, requires a stress-free lattice, and persists when anharmonic effects are included ([Aseginolaza et al., Nat. Phys. 20, 1288 (2024)](https://doi.org/10.1038/s41567-024-02441-z)). The default `asr="crystal"` imposes only translational invariance, so small numerical errors can leave the ZA branch linear or slightly imaginary near Gamma. For 2D materials:

- Relax the in-plane lattice with the same settings as the SCF, for example with `relax_job(..., relax_cell=True)` and `input_data={"cell": {"cell_dofree": "2Dxy"}}`.
- Use `asr="all"` in `matdyn.x`. It also imposes rotational invariance and the Huang conditions for vanishing stress (C. Lin, S. Poncé, N. Marzari, npj Comput. Mater. 8, 236 (2022)). For infrared-active layers, also set `write_lr` and `read_lr`, as described above.
- Check that `ω_ZA / |q|²` is nearly constant at small `|q|`.

```python
job_params = {
    "phonon_job": {"input_data": {"inputph": {"nq1": 6, "nq2": 6, "nq3": 1}}},
    "matdyn_job": {"input_data": {"input": {"asr": "all"}}},
}
```

For polar layers, `loto_2d=True` in both `q2r_job` and `matdyn_job` selects QE's two-dimensional treatment of the dipole-dipole interaction, in which the LO-TO splitting vanishes at Gamma.

## Reference data

Unmodified `matdyn.x` output for the Si example (QE 7.5, 3×3×3 DFPT grid) and a manifest of its inputs are stored under `tests/core/recipes/espresso_recipes/data/si_phonon_dispersion`, so that parsing and path bookkeeping are tested without a QE executable. The files were computed from the Latimer–Munro endpoints in band form, without the repeated interior Gamma points of the flow; Si has no LO-TO splitting, so the frequencies are the same. The files check the software, not the physics: the grid is not converged, and Si does not test the LO-TO handling.

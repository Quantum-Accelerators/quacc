"""FHI-aims recipes for the MAD-1.5 and MAD-1.6 datasets."""

from __future__ import annotations

from contextlib import contextmanager
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import TYPE_CHECKING

from ase.data import atomic_numbers

from quacc import Remove, get_settings, job
from quacc.recipes.aims._base import run_and_summarize

if TYPE_CHECKING:
    from collections.abc import Iterator
    from typing import Any

    from ase.atoms import Atoms

    from quacc.types import CopyFiles, RunSchema


@job
def static_job(
    atoms: Atoms,
    species_dir: str | Path = "tight",
    copy_files: CopyFiles | None = None,
    additional_fields: dict[str, Any] | None = None,
    **calc_kwargs: Any,
) -> RunSchema:
    """
    Run a single-point calculation with the shared MAD-1.5/1.6 r2SCAN protocol.

    Parameters
    ----------
    atoms
        The Atoms object. For the ``mc2d`` and ``mc3d_surface`` subsets, set
        ``atoms.info["subset"]`` accordingly to enable the dipole correction,
        or pass ``use_dipole_correction=True`` explicitly. Periodic structures
        must include a complete cell, with vacuum for slabs and 2D materials.
    species_dir
        Directory of FHI-aims 2020 species defaults, relative to
        ``AIMS_SPECIES_DEFAULTS`` or an absolute path. Defaults to ``"tight"``.
        The confined 5d function for Pr-Yb and 6d function for Pu-No are
        removed from temporary copies; the installed species files are unchanged.
    copy_files
        Files to copy (and decompress) from source to the runtime directory.
    additional_fields
        Additional fields to add to the results dictionary.
    **calc_kwargs
        Custom kwargs for the [ase.calculators.aims.Aims][] calculator.
        Set a value to ``quacc.Remove`` to remove a default parameter.

    Returns
    -------
    RunSchema
        Dictionary of results from [quacc.schemas.ase.Summarize.run][].

    Notes
    -----
    The published MAD-1.5 and MAD-1.6 ``control.in`` files are identical:
    https://archive.materialscloud.org/records/18tke-tt476 and
    https://archive.materialscloud.org/records/ffd5a-1hc37. The species
    modifications are described in the accompanying READMEs and
    https://arxiv.org/abs/2603.02089. The paper used FHI-aims version 250806.
    These recipes compute reference energies, forces, and periodic stresses;
    they do not generate structures or perform the dataset's outlier filtering.
    """
    periodic = bool(atoms.pbc.any())
    calc_defaults: dict[str, Any] = {
        "override_warning_libxc": True,
        "xc": "libxc MGGA_X_R2SCAN+MGGA_C_R2SCAN",
        "spin": "none",
        "charge": 0.0,
        "relativistic": "atomic_zora scalar",
        "sc_iter_limit": 200,
        "sc_init_iter": 99,
        "basis_threshold": 1e-5,
        "override_illconditioning": True,
        "override_kgrid_checks": True,
        "adjust_scf": "always 2",
        # ASE writes arbitrary string values verbatim. Both repeated keywords
        # are required to keep the 0.05 eV width in either adaptive SCF branch.
        "adjust_scf_param": (
            "lowgap occupation_width 0.05\n"
            "adjust_scf_param largegap occupation_width 0.05"
        ),
        "sc_accuracy_etot": 1e-6,
        "sc_accuracy_forces": 1e-4,
        "sc_accuracy_rho": 1e-5,
        "compute_forces": True,
    }
    if periodic:
        calc_defaults |= {"k_grid_density": 8, "compute_analytical_stress": True}
        if any(
            key in calc_kwargs and calc_kwargs[key] is not Remove
            for key in ("k_grid", "kpts", "kspacing")
        ):
            calc_defaults.pop("k_grid_density")
        if atoms.info.get("subset") in {"mc2d", "mc3d_surface"}:
            calc_defaults["use_dipole_correction"] = True

    with _species_directory(atoms, species_dir) as directory:
        calc_defaults["species_dir"] = directory
        return run_and_summarize(
            atoms,
            calc_defaults=calc_defaults,
            calc_swaps=calc_kwargs,
            additional_fields={"name": "MAD-1.5/1.6 Static"}
            | (additional_fields or {}),
            copy_files=copy_files,
        )


@contextmanager
def _species_directory(atoms: Atoms, species_dir: str | Path) -> Iterator[Path]:
    """Apply the MAD confined-orbital exclusions without editing installed files."""
    source = (get_settings().AIMS_SPECIES_DEFAULTS / species_dir).resolve()
    numbers = {atomic_numbers[symbol] for symbol in atoms.symbols}
    if not any(59 <= number <= 70 or 94 <= number <= 102 for number in numbers):
        yield source
        return

    with TemporaryDirectory(prefix="quacc-mad-") as tmpdir:
        directory = Path(tmpdir)
        for symbol in set(atoms.symbols):
            number = atomic_numbers[symbol]
            filename = f"{number:02}_{symbol}_default"
            text = (source / filename).read_text()
            shell = "5" if 59 <= number <= 70 else "6" if 94 <= number <= 102 else None
            if shell:
                text = "".join(
                    line
                    for line in text.splitlines(keepends=True)
                    if line.split()[:3] != ["confined", shell, "d"]
                )
            (directory / filename).write_text(text)
        yield directory

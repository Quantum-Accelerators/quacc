"""Base jobs for FHI-aims."""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

from ase.calculators.aims import Aims, AimsProfile

from quacc import get_settings
from quacc.runners.ase import Runner
from quacc.schemas.ase import Summarize
from quacc.utils.dicts import recursive_dict_merge
from quacc.utils.kpts import kspacing_to_grid

if TYPE_CHECKING:
    from typing import Any

    from ase.atoms import Atoms

    from quacc.types import Filenames, OptParams, RunSchema, SourceDirectory


def run_and_summarize(
    atoms: Atoms,
    calc_defaults: dict[str, Any] | None = None,
    calc_swaps: dict[str, Any] | None = None,
    additional_fields: dict[str, Any] | None = None,
    copy_files: SourceDirectory | dict[SourceDirectory, Filenames] | None = None,
) -> RunSchema:
    """
    Base function to carry out FHI-aims recipes.

    Parameters
    ----------
    atoms
        Atoms object
    calc_defaults
        The default calculator parameters.
    calc_swaps
        Custom kwargs for the FHI-aims calculator. Set a value to
        `quacc.Remove` to remove a pre-existing key entirely. For a list of available
        keys, refer to the [ase.calculators.aims.Aims][] calculator.
    additional_fields
        Any additional fields to supply to the summarizer.
    copy_files
        Files to copy (and decompress) from source to the runtime directory.

    Returns
    -------
    RunSchema
        Dictionary of results from [quacc.schemas.ase.Summarize.run][]
    """
    calc = prep_calculator(atoms, calc_defaults=calc_defaults, calc_swaps=calc_swaps)
    final_atoms = Runner(atoms, calc, copy_files=copy_files).run_calc()

    return Summarize(move_magmoms=True, additional_fields=additional_fields).run(
        final_atoms, atoms
    )


def run_and_summarize_opt(
    atoms: Atoms,
    calc_defaults: dict[str, Any] | None = None,
    calc_swaps: dict[str, Any] | None = None,
    opt_defaults: dict[str, Any] | None = None,
    opt_params: OptParams | None = None,
    additional_fields: dict[str, Any] | None = None,
    copy_files: SourceDirectory | dict[SourceDirectory, Filenames] | None = None,
) -> RunSchema:
    """
    Base function to carry out FHI-aims recipes with ASE optimizers.

    Parameters
    ----------
    atoms
        Atoms object
    calc_defaults
        The default calculator parameters.
    calc_swaps
        Custom kwargs for the FHI-aims calculator. Set a value to
        `quacc.Remove` to remove a pre-existing key entirely. For a list of available
        keys, refer to the [ase.calculators.aims.Aims][] calculator.
    opt_defaults
        The default optimization parameters.
    opt_params
        Dictionary of custom kwargs for the optimization process. For a list
        of available keys, refer to [quacc.runners.ase.Runner.run_opt][].
    additional_fields
        Any additional fields to supply to the summarizer.
    copy_files
        Files to copy (and decompress) from source to the runtime directory.

    Returns
    -------
    RunSchema
        Dictionary of results from [quacc.schemas.ase.Summarize.run][]
    """
    opt_flags = recursive_dict_merge(opt_defaults, opt_params)
    calc = prep_calculator(atoms, calc_defaults=calc_defaults, calc_swaps=calc_swaps)
    dyn = Runner(atoms, calc, copy_files=copy_files).run_opt(**opt_flags)

    return Summarize(move_magmoms=True, additional_fields=additional_fields).opt(dyn)


def prep_calculator(
    atoms: Atoms,
    calc_defaults: dict[str, Any] | None = None,
    calc_swaps: dict[str, Any] | None = None,
) -> Aims:
    """
    Prepare the FHI-aims calculator.

    Parameters
    ----------
    atoms
        Atoms object
    calc_defaults
        The default calculator parameters.
    calc_swaps
        Custom kwargs for the FHI-aims calculator. Set a value to
        `quacc.Remove` to remove a pre-existing key entirely. For a list of available
        keys, refer to the [ase.calculators.aims.Aims][] calculator.

    Returns
    -------
    Aims
        The FHI-aims calculator.
    """
    calc_flags = recursive_dict_merge(calc_defaults or {}, calc_swaps or {})
    settings = get_settings()
    species_dir = calc_flags.pop("species_dir", None)

    if not any(atoms.pbc):
        for key in ["kspacing", "k_grid", "k_grid_density"]:
            if key in calc_flags:
                calc_flags.pop(key)
    elif (
        "kspacing" in calc_flags
        and "k_grid" not in calc_flags
        and "k_grid_density" not in calc_flags
    ):
        kspacing = calc_flags.pop("kspacing")
        calc_flags["k_grid"] = kspacing_to_grid(atoms, kspacing)

    if "spin" not in calc_flags and hasattr(atoms, "get_initial_magnetic_moments"):
        magmoms = atoms.get_initial_magnetic_moments()
        if magmoms is not None and any(abs(m) > 1e-6 for m in magmoms):
            calc_flags["spin"] = "collinear"

    aims_cmd = f"{settings.AIMS_PARALLEL_CMD} {settings.AIMS_BIN}"
    species_path = settings.AIMS_SPECIES_DEFAULTS
    if species_dir:
        species_path = Path(species_path) / species_dir

    return Aims(
        profile=AimsProfile(
            command=aims_cmd.strip(), default_species_directory=str(species_path)
        ),
        **calc_flags,
    )

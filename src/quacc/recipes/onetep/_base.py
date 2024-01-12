"""Base jobs for Onetep."""
from __future__ import annotations

from typing import TYPE_CHECKING

from ase import Atoms

### Do we need to import from quacc at this point we can do it directly from ase
from quacc.calculators.onetep.onetep import Onetep, OnetepProfile
from ase.calculators.onetep import OnetepTemplate
from quacc.runners.ase import run_calc
from quacc.schemas.ase import summarize_run

if TYPE_CHECKING:
    from typing import Any

    from quacc.schemas._aliases.ase import RunSchema


def base_fn(
    atoms: Atoms = Atoms(),
    preset: str | None = None,
    template: OnetepTemplate | None = None,
    profile: OnetepProfile | None = None,
    calc_defaults: dict[str, Any] | None = None,
    calc_swaps: dict[str, Any] | None = None,
    parallel_info: dict[str] | None = None,
    additional_fields: dict[str, Any] | None = None,
    copy_files: list[str] | None = None,
) -> RunSchema:
    """
    Base function to carry out Onetep recipes.

    Parameters
    ----------
    atoms
        Atoms object
    preset
        Name of the preset to use
    template
        OnetepTemplate to use
    profile
        OnetepProfile to use
    calc_defaults
        The default calculator parameters.
    calc_swaps
        Custom kwargs for the Onetep calculator. Set a value to
        `None` to remove a pre-existing key entirely. For a list of available
        keys, refer to the `ase.calculators.Onetep.Onetep` calculator.
    parallel_info
        Dictionary of parallelization information.
    additional_fields
        Any additional fields to supply to the summarizer.
    copy_files
        Files to copy to the runtime directory.

    Returns
    -------
    RunSchema
        Dictionary of results from [quacc.schemas.ase.summarize_run][]
    """

    atoms.calc = Onetep(
        input_atoms=atoms,
        preset=preset,
        template=template,
        profile=profile,
        calc_defaults=calc_defaults,
        parallel_info=parallel_info,
        **calc_swaps,
    )

    final_atoms = run_calc(atoms, copy_files=copy_files)

    return summarize_run(
        final_atoms, input_atoms=atoms, additional_fields=additional_fields
    )

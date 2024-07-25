"""Base jobs for MRCC."""

from __future__ import annotations

from typing import TYPE_CHECKING

from quacc.calculators.mrcc import MRCC
from quacc.runners.ase import Runner
from quacc.schemas.ase import summarize_run
from quacc.utils.dicts import recursive_dict_merge

if TYPE_CHECKING:
    from typing import Any

    from ase.atoms import Atoms

    from quacc.types import Filenames, OptParams, OptSchema, RunSchema, SourceDirectory


def run_and_summarize(
    atoms: Atoms,
    charge: int = 0,
    spin_multiplicity: int = 1,
    calc_defaults: dict[str, Any] | None = None,
    calc_swaps: dict[str, Any] | None = None,
    additional_fields: dict[str, Any] | None = None,
    copy_files: SourceDirectory | dict[SourceDirectory, Filenames] | None = None,
) -> RunSchema:
    """
    Base job function used for MRCC recipes that don't rely on ASE optimizers or other
    ASE dynamics classes.

    Parameters
    ----------
    atoms
        Atoms object
    charge
        Charge of the system.
    spin_multiplicity
        Multiplicity of the system.
    calc_defaults
        The default parameters for the recipe.
    calc_swaps
        Dictionary of custom kwargs for the MRCC calculator. Set a value to `quacc.Remove` to
        remove a pre-existing key entirely. For a list of available keys, refer to the
        `quacc.calculators.mrcc.mrcc.MRCC` calculator.
    additional_fields
        Any additional fields to set in the summary.
    copy_files
        Files to copy (and decompress) from source to the runtime directory.

    Returns
    -------
    RunSchema
        Dictionary of results from [quacc.schemas.ase.summarize_run][]
    """
    calc_flags = recursive_dict_merge(calc_defaults, calc_swaps)
    calc = MRCC(
        atoms, charge=charge, spin_multiplicity=spin_multiplicity, **calc_flags
    )
    final_atoms = Runner(atoms, calc, copy_files=copy_files).run_calc()

    return summarize_run(
        final_atoms,
        atoms,
        charge_and_multiplicity=(charge, spin_multiplicity),
        additional_fields=additional_fields,
    )
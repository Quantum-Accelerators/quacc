"""Base jobs for Psi4."""

from __future__ import annotations

from typing import TYPE_CHECKING

from ase.calculators.psi4 import Psi4

from quacc.runners.ase import Runner
from quacc.schemas.ase import Summarize
from quacc.utils.dicts import recursive_dict_merge

if TYPE_CHECKING:
    from typing import Any

    from ase.atoms import Atoms

    from quacc.types import Filenames, RunSchema, SourceDirectory


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
    Base function to carry out Psi4 recipes.

    Parameters
    ----------
    atoms
        Atoms object
    charge
        Charge of the system.
    spin_multiplicity
        Multiplicity of the system.
    calc_defaults
        The default calculator parameters.
    calc_swaps
        Custom kwargs for the Psi4 calculator. Set a value to
        `quacc.Remove` to remove a pre-existing key entirely. For a list of available
        keys, refer to the [ase.calculators.psi4.Psi4][] calculator.
    additional_fields
        Any additional fields to supply to the summarizer.
    copy_files
        Files to copy (and decompress) from source to the runtime directory.

    Returns
    -------
    RunSchema
        Dictionary of results from [quacc.schemas.ase.Summarize.run][]
    """
    calc_flags = recursive_dict_merge(calc_defaults, calc_swaps)

    calc = Psi4(**calc_flags)
    final_atoms = Runner(atoms, calc, copy_files=copy_files).run_calc()

    return Summarize(
        charge_and_multiplicity=(charge, spin_multiplicity),
        additional_fields=additional_fields,
    ).run(final_atoms, atoms)

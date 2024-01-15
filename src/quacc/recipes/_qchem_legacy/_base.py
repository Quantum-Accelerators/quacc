"""Base jobs for Q-Chem."""
from __future__ import annotations

from typing import TYPE_CHECKING

from quacc.calculators._qchem_legacy import QChem
from quacc.runners.ase import run_calc, run_opt
from quacc.schemas.ase import summarize_opt_run, summarize_run
from quacc.utils.dicts import recursive_dict_merge, remove_dict_entries

if TYPE_CHECKING:
    from pathlib import Path
    from typing import Any

    from ase.atoms import Atoms

    from quacc.schemas._aliases.ase import OptSchema, RunSchema


def base_fn(
    atoms: Atoms,
    charge: int = 0,
    spin_multiplicity: int = 1,
    defaults: dict[str, Any] | None = None,
    additional_fields: dict[str, Any] | None = None,
    copy_files: str | Path | list[str | Path] | None = None,
) -> RunSchema:
    """
    Base job function used for Q-Chem recipes that don't rely on ASE optimizers or other
    ASE dynamics classes.

    Parameters
    ----------
    atoms
        Atoms object
    charge
        Charge of the system.
    spin_multiplicity
        Multiplicity of the system.
    defaults
        The default parameters for the recipe.
    additional_fields
        Any additional fields to set in the summary.
    copy_files
        File(s) to copy to the runtime directory. If a directory is provided, it will be recursively unpacked.

    Returns
    -------
    RunSchema
        Dictionary of results from [quacc.schemas.ase.summarize_run][]
    """

    qchem_flags = remove_dict_entries(defaults, None)

    atoms.calc = QChem(atoms, **qchem_flags)
    final_atoms = run_calc(atoms, copy_files=copy_files)

    return summarize_run(
        final_atoms,
        atoms,
        charge_and_multiplicity=(charge, spin_multiplicity),
        additional_fields=additional_fields,
    )


def base_opt_fn(
    atoms: Atoms,
    charge: int = 0,
    spin_multiplicity: int = 1,
    calc_defaults: dict[str, Any] | None = None,
    opt_defaults: dict[str, Any] | None = None,
    opt_params: dict[str, Any] | None = None,
    additional_fields: dict[str, Any] | None = None,
    copy_files: str | Path | list[str | Path] | None = None,
) -> OptSchema:
    """
    Base function for Q-Chem recipes that involve ASE optimizers.

    Parameters
    ----------
    atoms
        Atoms object
    charge
        Charge of the system.
    spin_multiplicity
        Multiplicity of the system.
    calc_defaults
        Default arguments for the Q-Chem calculator.
    opt_defaults
        Default arguments for the ASE optimizer.
    opt_params
        Dictionary of custom kwargs for [quacc.runners.ase.run_opt][]
    copy_files
        File(s) to copy to the runtime directory. If a directory is provided, it will be recursively unpacked.

    Returns
    -------
    OptSchema
        Dictionary of results from [quacc.schemas.ase.summarize_opt_run][]
    """
    # TODO:
    #   - passing initial Hessian?

    qchem_flags = remove_dict_entries(calc_defaults, None)
    opt_flags = recursive_dict_merge(opt_defaults, opt_params)

    atoms.calc = QChem(atoms, **qchem_flags)
    dyn = run_opt(atoms, copy_files=copy_files, **opt_flags)

    return summarize_opt_run(
        dyn,
        charge_and_multiplicity=(charge, spin_multiplicity),
        additional_fields=additional_fields,
    )

"""Base jobs for ORCA."""

from __future__ import annotations

from typing import TYPE_CHECKING

from ase.calculators.orca import ORCA, OrcaProfile, OrcaTemplate

from quacc import SETTINGS
from quacc.atoms.core import get_final_atoms_from_dyn
from quacc.runners.ase import run_calc, run_opt
from quacc.schemas.ase import summarize_opt_run
from quacc.schemas.cclib import cclib_summarize_run
from quacc.utils.dicts import recursive_dict_merge
from quacc.utils.lists import merge_list_params

if TYPE_CHECKING:
    from typing import Any

    from ase.atoms import Atoms

    from quacc.schemas._aliases.cclib import cclibASEOptSchema, cclibSchema
    from quacc.utils.files import Filenames, SourceDirectory

_LABEL = OrcaTemplate()._label  # skipcq: PYL-W0212
LOG_FILE = f"{_LABEL}.out"
GEOM_FILE = f"{_LABEL}.xyz"


def base_fn(
    atoms: Atoms,
    charge: int = 0,
    spin_multiplicity: int = 1,
    default_inputs: list[str] | None = None,
    default_blocks: list[str] | None = None,
    input_swaps: list[str] | None = None,
    block_swaps: list[str] | None = None,
    additional_fields: dict[str, Any] | None = None,
    copy_files: SourceDirectory | dict[SourceDirectory, Filenames] | None = None,
) -> cclibSchema:
    """
    Base job function for ORCA recipes.

    Parameters
    ----------
    atoms
        Atoms object
    charge
        Charge of the system.
    spin_multiplicity
        Multiplicity of the system.
    default_inputs
        Default input parameters.
    default_blocks
        Default block input parameters.
    input_swaps
        List of orcasimpleinput swaps for the calculator. To remove entries
        from the defaults, put a `#` in front of the name.
    block_swaps
        List of orcablock swaps for the calculator. To remove entries
        from the defaults, put a `#` in front of the name.
    additional_fields
        Any additional fields to supply to the summarizer.
    copy_files
        Files to copy (and decompress) from source to the runtime directory.

    Returns
    -------
    cclibSchema
        Dictionary of results from [quacc.schemas.cclib.cclib_summarize_run][]
    """
    inputs = merge_list_params(default_inputs, input_swaps)
    blocks = merge_list_params(default_blocks, block_swaps)
    orcasimpleinput = " ".join(inputs)
    orcablocks = " ".join(blocks)

    atoms.calc = ORCA(
        profile=OrcaProfile(SETTINGS.ORCA_CMD),
        charge=charge,
        mult=spin_multiplicity,
        orcasimpleinput=orcasimpleinput,
        orcablocks=orcablocks,
    )
    atoms = run_calc(atoms, geom_file=GEOM_FILE, copy_files=copy_files)

    return cclib_summarize_run(atoms, LOG_FILE, additional_fields=additional_fields)


def base_opt_fn(
    atoms: Atoms,
    charge: int = 0,
    spin_multiplicity: int = 1,
    default_inputs: list[str] | None = None,
    default_blocks: list[str] | None = None,
    input_swaps: list[str] | None = None,
    block_swaps: list[str] | None = None,
    opt_defaults: dict[str, Any] | None = None,
    opt_params: dict[str, Any] | None = None,
    additional_fields: dict[str, Any] | None = None,
    copy_files: SourceDirectory | dict[SourceDirectory, Filenames] | None = None,
) -> cclibASEOptSchema:
    """
    Base job function for ORCA recipes with ASE optimizer.

    Parameters
    ----------
    atoms
        Atoms object
    charge
        Charge of the system.
    spin_multiplicity
        Multiplicity of the system.
    default_inputs
        Default input parameters.
    default_blocks
        Default block input parameters.
    input_swaps
        List of orcasimpleinput swaps for the calculator. To remove entries
        from the defaults, put a `#` in front of the name.
    block_swaps
        List of orcablock swaps for the calculator. To remove entries
        from the defaults, put a `#` in front of the name.
    opt_defaults
        Default arguments for the ASE optimizer.
    opt_params
        Dictionary of custom kwargs for [quacc.runners.ase.run_opt][]
    additional_fields
        Any additional fields to supply to the summarizer.
    copy_files
        Files to copy (and decompress) from source to the runtime directory.

    Returns
    -------
    cclibASEOptSchema
        Dictionary of results from [quacc.schemas.cclib.cclib_summarize_run][]
    """
    inputs = merge_list_params(default_inputs, input_swaps)
    blocks = merge_list_params(default_blocks, block_swaps)
    opt_flags = recursive_dict_merge(opt_defaults, opt_params)
    orcasimpleinput = " ".join(inputs)
    orcablocks = " ".join(blocks)

    atoms.calc = ORCA(
        profile=OrcaProfile(SETTINGS.ORCA_CMD),
        charge=charge,
        mult=spin_multiplicity,
        orcasimpleinput=orcasimpleinput,
        orcablocks=orcablocks,
    )

    dyn = run_opt(atoms, copy_files=copy_files, **opt_flags)
    opt_run_summary = summarize_opt_run(
        dyn,
        charge_and_multiplicity=(charge, spin_multiplicity),
        additional_fields=additional_fields,
    )

    final_atoms = get_final_atoms_from_dyn(dyn)
    cclib_summary = cclib_summarize_run(
        final_atoms, LOG_FILE, additional_fields=additional_fields
    )
    return recursive_dict_merge(cclib_summary, opt_run_summary)

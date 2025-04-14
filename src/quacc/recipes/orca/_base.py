"""Base jobs for ORCA."""

from __future__ import annotations

from typing import TYPE_CHECKING

from ase.calculators.orca import ORCA, OrcaProfile, OrcaTemplate

from quacc import get_settings
from quacc.runners.ase import Runner
from quacc.schemas.ase import Summarize
from quacc.utils.dicts import recursive_dict_merge
from quacc.utils.lists import merge_list_params

if TYPE_CHECKING:
    from typing import Any

    from ase.atoms import Atoms

    from quacc.types import Filenames, OptParams, OptSchema, RunSchema, SourceDirectory

_LABEL = OrcaTemplate()._label  # skipcq: PYL-W0212
GEOM_FILE = f"{_LABEL}.xyz"


def run_and_summarize(
    atoms: Atoms,
    charge: int = 0,
    spin_multiplicity: int = 1,
    default_inputs: list[str] | None = None,
    default_blocks: list[str] | None = None,
    input_swaps: list[str] | None = None,
    block_swaps: list[str] | None = None,
    additional_fields: dict[str, Any] | None = None,
    copy_files: SourceDirectory | dict[SourceDirectory, Filenames] | None = None,
    **calc_kwargs,
) -> RunSchema:
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
    **calc_kwargs
        Any other keyword arguments to pass to the `ORCA` calculator.

    Returns
    -------
    RunSchema
        Dictionary of results
    """
    calc = prep_calculator(
        charge=charge,
        spin_multiplicity=spin_multiplicity,
        default_inputs=default_inputs,
        default_blocks=default_blocks,
        input_swaps=input_swaps,
        block_swaps=block_swaps,
        **calc_kwargs,
    )

    final_atoms = Runner(atoms, calc, copy_files=copy_files).run_calc(
        geom_file=GEOM_FILE
    )

    return Summarize(additional_fields=additional_fields).run(final_atoms, atoms)


def run_and_summarize_opt(
    atoms: Atoms,
    charge: int = 0,
    spin_multiplicity: int = 1,
    default_inputs: list[str] | None = None,
    default_blocks: list[str] | None = None,
    input_swaps: list[str] | None = None,
    block_swaps: list[str] | None = None,
    opt_defaults: dict[str, Any] | None = None,
    opt_params: OptParams | None = None,
    additional_fields: dict[str, Any] | None = None,
    copy_files: SourceDirectory | dict[SourceDirectory, Filenames] | None = None,
    **calc_kwargs,
) -> OptSchema:
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
        Dictionary of custom kwargs for [quacc.runners.ase.Runner.run_opt][]
    additional_fields
        Any additional fields to supply to the summarizer.
    copy_files
        Files to copy (and decompress) from source to the runtime directory.
    **calc_kwargs
        Any other keyword arguments to pass to the `ORCA` calculator.

    Returns
    -------
    OptSchema
        Dictionary of results
    """
    calc = prep_calculator(
        charge=charge,
        spin_multiplicity=spin_multiplicity,
        default_inputs=default_inputs,
        default_blocks=default_blocks,
        input_swaps=input_swaps,
        block_swaps=block_swaps,
        **calc_kwargs,
    )

    opt_flags = recursive_dict_merge(opt_defaults, opt_params)
    dyn = Runner(atoms, calc, copy_files=copy_files).run_opt(**opt_flags)
    return Summarize(additional_fields=additional_fields).opt(dyn)


def prep_calculator(
    charge: int = 0,
    spin_multiplicity: int = 1,
    default_inputs: list[str] | None = None,
    default_blocks: list[str] | None = None,
    input_swaps: list[str] | None = None,
    block_swaps: list[str] | None = None,
    **calc_kwargs,
) -> ORCA:
    """
    Prepare the ORCA calculator.

    Parameters
    ----------
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
    **calc_kwargs
        Any other keyword arguments to pass to the `ORCA` calculator.

    Returns
    -------
    ORCA
        The ORCA calculator
    """
    inputs = merge_list_params(default_inputs, input_swaps)
    blocks = merge_list_params(default_blocks, block_swaps)
    if "xyzfile" not in inputs:
        inputs.append("xyzfile")
    orcasimpleinput = " ".join(inputs)
    orcablocks = "\n".join(blocks)
    settings = get_settings()

    return ORCA(
        profile=OrcaProfile(command=settings.ORCA_CMD),
        charge=charge,
        mult=spin_multiplicity,
        orcasimpleinput=orcasimpleinput,
        orcablocks=orcablocks,
        **calc_kwargs,
    )

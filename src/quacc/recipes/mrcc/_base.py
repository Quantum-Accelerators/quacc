"""Base jobs for MRCC."""

from __future__ import annotations

from typing import TYPE_CHECKING

from quacc import get_settings
from quacc.calculators.mrcc.mrcc import MRCC, MrccProfile
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
    default_inputs: dict[str, str] | None = None,
    blocks: str | None = None,
    input_swaps: dict[str, str] | None = None,
    additional_fields: dict[str, Any] | None = None,
    copy_files: SourceDirectory | dict[SourceDirectory, Filenames] | None = None,
    **calc_kwargs,
) -> RunSchema:
    """
    Base job function for MRCC recipes.

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
    blocks
        Block input parameters.
    input_swaps
        List of orcasimpleinput swaps for the calculator. To remove entries
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
        Dictionary of results from [quacc.schemas.ase.Summarize][]
    """
    calc = prep_calculator(
        charge=charge,
        spin_multiplicity=spin_multiplicity,
        default_inputs=default_inputs,
        blocks=blocks,
        input_swaps=input_swaps,
        **calc_kwargs,
    )

    final_atoms = Runner(atoms, calc, copy_files=copy_files).run_calc()

    return Summarize(
        charge_and_multiplicity=(charge, spin_multiplicity),
        additional_fields=additional_fields,
    ).run(final_atoms, atoms)


def prep_calculator(
    charge: int = 0,
    spin_multiplicity: int = 1,
    default_inputs: dict[str, str] | None = None,
    blocks: str | None = None,
    input_swaps: dict[str, str] | None = None,
    **calc_kwargs,
) -> MRCC:
    """
    Prepare the MRCC calculator.

    Parameters
    ----------
    charge
        Charge of the system.
    spin_multiplicity
        Multiplicity of the system.
    default_inputs
        Default input parameters.
    blocks
        MRCC block input string.
    input_swaps
        List of mrccinput swaps for the calculator. To remove entries
        from the defaults, put a `#` in front of the name.
    **calc_kwargs
        Any other keyword arguments to pass to the `ORCA` calculator.

    Returns
    -------
    MRCC
        The MRCC calculator
    """
    mrccinput = recursive_dict_merge(default_inputs, input_swaps)

    # Check that there is no duplication of keywords between mrccinput and blocks
    if blocks is not None:
        for line in blocks.split():
            if "=" in line:
                key = line.split("=")[0]
                if key == "scftype":
                    raise ValueError(
                        "Keyword scftype must be specified in mrccinput only"
                    )
                if key in mrccinput:
                    raise ValueError(
                        f"Keyword {key} is duplicated in both mrccinput and blocks"
                    )

    # If spin_multiplicity bigger than 1, check if scftype is in either mrccinput or blocks
    if spin_multiplicity > 1 and "scftype" not in mrccinput:
        raise ValueError(
            "For spin_multiplicity > 1, scftype keyword must be specified in mrccinput"
        )
    elif spin_multiplicity > 1 and mrccinput["scftype"].lower() not in [
        "uhf",
        "uks",
        "rohf",
        "roks",
    ]:
        raise ValueError(
            "For spin_multiplicity > 1, scftype must not be set to RHF or RKS"
        )

    settings = get_settings()

    return MRCC(
        profile=MrccProfile(command=settings.MRCC_CMD),
        charge=charge,
        mult=spin_multiplicity,
        mrccinput=mrccinput,
        mrccblocks=blocks,
        **calc_kwargs,
    )

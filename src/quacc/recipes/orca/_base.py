"""Base jobs for ORCA"""
from __future__ import annotations

from typing import TYPE_CHECKING

from ase.calculators.orca import ORCA, OrcaProfile, OrcaTemplate

from quacc import SETTINGS
from quacc.runners.ase import run_calc
from quacc.schemas.cclib import cclib_summarize_run
from quacc.utils.dicts import recursive_dict_merge

if TYPE_CHECKING:
    from pathlib import Path
    from typing import Any

    from ase.atoms import Atoms

    from quacc.schemas._aliases.cclib import cclibSchema

_LABEL = OrcaTemplate()._label  # skipcq: PYL-W0212
LOG_FILE = f"{_LABEL}.out"
GEOM_FILE = f"{_LABEL}.xyz"


def base_fn(
    atoms: Atoms,
    charge: int = 0,
    spin_multiplicity: int = 1,
    default_inputs: dict[str, Any] | None = None,
    default_blocks: dict[str, Any] | None = None,
    input_swaps: dict[str, Any] | None = None,
    block_swaps: dict[str, Any] | None = None,
    additional_fields: dict[str, Any] | None = None,
    copy_files: str | Path | list[str | Path] | None = None,
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
        Dictionary of orcasimpleinput swaps for the calculator. To enable new
        entries, set the value as True. To remove entries from the defaults, set
        the value as None.
    block_swaps
        Dictionary of orcablock swaps for the calculator. To enable new entries,
        set the value as True. To remove entries from the defaults, set the
        value as None.
    additional_fields
        Any additional fields to supply to the summarizer.
    copy_files
        File(s) to copy to the runtime directory. If a directory is provided, it will be recursively unpacked.

    Returns
    -------
    cclibSchema
        Dictionary of results from [quacc.schemas.cclib.cclib_summarize_run][]
    """
    inputs = recursive_dict_merge(default_inputs, input_swaps)
    blocks = recursive_dict_merge(default_blocks, block_swaps)
    orcasimpleinput = " ".join(list(inputs.keys()))
    orcablocks = " ".join(list(blocks.keys()))

    atoms.calc = ORCA(
        profile=OrcaProfile(SETTINGS.ORCA_CMD),
        charge=charge,
        mult=spin_multiplicity,
        orcasimpleinput=orcasimpleinput,
        orcablocks=orcablocks,
    )
    atoms = run_calc(atoms, geom_file=GEOM_FILE, copy_files=copy_files)

    return cclib_summarize_run(atoms, LOG_FILE, additional_fields=additional_fields)

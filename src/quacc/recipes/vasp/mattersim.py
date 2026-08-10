"""MatterSim-compatible VASP recipes."""

from __future__ import annotations

from typing import TYPE_CHECKING

from pymatgen.io.vasp.sets import MPRelaxSet

from quacc import job
from quacc.calculators.vasp.params import MPtoASEConverter
from quacc.recipes.vasp._base import run_and_summarize

if TYPE_CHECKING:
    from typing import Any

    from ase.atoms import Atoms

    from quacc.types import SourceDirectory, VaspSchema
    from quacc.wflow_tools.job_argument import Copy


@job
def mattersim_static_job(
    atoms: Atoms,
    copy_files: SourceDirectory | Copy | None = None,
    additional_fields: dict[str, Any] | None = None,
    **calc_kwargs,
) -> VaspSchema:
    """
    Run a static calculation with the VASP settings used for MatterSim.

    Parameters
    ----------
    atoms
        Atoms object.
    copy_files
        Files to copy (and decompress) from source to the runtime directory.
    additional_fields
        Additional fields to add to the results dictionary.
    **calc_kwargs
        Custom kwargs for the Vasp calculator. Set a value to `None` to remove a
        pre-existing key entirely. All ASE Vasp calculator keyword arguments are
        supported.

    Returns
    -------
    VaspSchema
        Dictionary of results from [quacc.schemas.vasp.VaspSummarize.run][].
    """
    calc_defaults = MPtoASEConverter(atoms=atoms).convert_input_set(MPRelaxSet())
    calc_defaults["nsw"] = 0

    return run_and_summarize(
        atoms,
        calc_defaults=calc_defaults,
        calc_swaps=calc_kwargs,
        additional_fields={"name": "MatterSim Static"} | (additional_fields or {}),
        copy_files=copy_files,
    )

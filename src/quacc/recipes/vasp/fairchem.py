"""
Meta FAIR recipes

!!! Important

    Calculations will be use the following version of pseudopotentials:
    - OMat: v.54
    - OMC: v.54
"""

from __future__ import annotations

from importlib.util import find_spec
from pathlib import Path
from typing import TYPE_CHECKING

from monty.dev import requires

from quacc import job
from quacc.calculators.vasp.params import MPtoASEConverter
from quacc.recipes.vasp._base import run_and_summarize

has_fairchem = bool(find_spec("fairchem"))
has_fairchem_omat = has_fairchem and bool(find_spec("fairchem.data.omat"))
has_fairchem_omc = has_fairchem and bool(find_spec("fairchem.data.omc"))
has_atomate2 = bool(find_spec("atomate2"))

if TYPE_CHECKING:
    from typing import Any

    from ase.atoms import Atoms

    from quacc.types import SourceDirectory, VaspSchema
    from quacc.wflow_tools.job_argument import Copy


@job
@requires(
    has_fairchem_omat,
    "fairchem-data-omat is not installed. Run `pip install quacc[fairchem]`",
)
def omat_static_job(
    atoms: Atoms,
    copy_files: SourceDirectory | Copy | None = None,
    additional_fields: dict[str, Any] | None = None,
    **calc_kwargs,
) -> VaspSchema:
    """
    Carry out a static calculation with OMat settings.

    Parameters
    ----------
    atoms
        Atoms object
    copy_files
        Files to copy (and decompress) from source to the runtime directory.
    additional_fields
        Additional fields to add to the results dictionary.
    **calc_kwargs
        Custom kwargs for the Vasp calculator. Set a value to
        `None` to remove a pre-existing key entirely. For a list of available
        keys, refer to [quacc.calculators.vasp.vasp.Vasp][]. All of the ASE
        Vasp calculator keyword arguments are supported.

    Returns
    -------
    VaspSchema
        Dictionary of results from [quacc.schemas.vasp.VaspSummarize.run][].
        See the type-hint for the data structure.
    """
    from fairchem.data.omat.vasp.sets import OMat24StaticSet

    calc_defaults = MPtoASEConverter(atoms=atoms).convert_input_set(
        OMat24StaticSet()
    ) | {"incar_copilot": "ncore"}

    return run_and_summarize(
        atoms,
        calc_defaults=calc_defaults,
        calc_swaps=calc_kwargs,
        additional_fields={"name": "OMat Static"} | (additional_fields or {}),
        copy_files=copy_files,
    )


@job
@requires(
    has_fairchem_omc,
    "fairchem-data-omc is not installed. Run `pip install quacc[fairchem]`",
)
@requires(has_atomate2, "atomate2 is not installed. Run `pip install quacc[fairchem]`")
def omc_static_job(
    atoms: Atoms,
    copy_files: SourceDirectory | Copy | None = None,
    additional_fields: dict[str, Any] | None = None,
    **calc_kwargs,
) -> VaspSchema:
    """
    Carry out a static calculation with OMC settings.

    Parameters
    ----------
    atoms
        Atoms object
    copy_files
        Files to copy (and decompress) from source to the runtime directory.
    additional_fields
        Additional fields to add to the results dictionary.
    **calc_kwargs
        Custom kwargs for the Vasp calculator. Set a value to
        `None` to remove a pre-existing key entirely. For a list of available
        keys, refer to [quacc.calculators.vasp.vasp.Vasp][]. All of the ASE
        Vasp calculator keyword arguments are supported.

    Returns
    -------
    VaspSchema
        Dictionary of results from [quacc.schemas.vasp.VaspSummarize.run][].
        See the type-hint for the data structure.
    """
    from argparse import Namespace

    from fairchem.data.omc.scripts import create_vasp_inputs

    incar_path = Namespace(
        incar_yml_dir=Path(create_vasp_inputs.__file__).parent / "incars", type="static"
    )
    input_generator = create_vasp_inputs.create_input_generator(incar_path)
    calc_defaults = MPtoASEConverter(atoms=atoms).convert_input_generator(
        input_generator
    ) | {"incar_copilot": "light"}

    return run_and_summarize(
        atoms,
        calc_defaults=calc_defaults,
        calc_swaps=calc_kwargs,
        additional_fields={"name": "OMC Static"} | (additional_fields or {}),
        copy_files=copy_files,
    )

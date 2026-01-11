"""
Meta FAIR recipes

!!! Important

    Make sure that you use the appropriate version pseudpotentials:
    - OMat: v.54
    - OMC: v.54
"""

from __future__ import annotations

from importlib.util import find_spec
from typing import TYPE_CHECKING

from monty.dev import requires

from quacc import job
from quacc.calculators.vasp.params import MPtoASEConverter
from quacc.recipes.vasp._base import run_and_summarize

has_fairchem = bool(find_spec("fairchem"))
has_fairchem_omat = has_fairchem and bool(find_spec("fairchem.data.omat"))
has_atomate2 = bool(find_spec("atomate2"))

if TYPE_CHECKING:
    from typing import Any

    from ase.atoms import Atoms

    from quacc.types import Filenames, SourceDirectory, VaspSchema


@job
@requires(
    has_fairchem_omat,
    "fairchem-data-omat is not installed. Run `pip install quacc[fairchem]`",
)
def omat_static_job(
    atoms: Atoms,
    copy_files: SourceDirectory | dict[SourceDirectory, Filenames] | None = None,
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

    calc_defaults = MPtoASEConverter(atoms=atoms).convert_input_set(OMat24StaticSet())
    calc_defaults |= {"pp_version": "54", "incar_copilot": False}

    return run_and_summarize(
        atoms,
        calc_defaults=calc_defaults,
        calc_swaps=calc_kwargs,
        additional_fields={"name": "OMat Static"} | (additional_fields or {}),
        copy_files=copy_files,
    )


@job
@requires(has_atomate2, "atomate2 is not installed. Run `pip install quacc[fairchem]`")
def omc_static_job(
    atoms: Atoms,
    copy_files: SourceDirectory | dict[SourceDirectory, Filenames] | None = None,
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

    calc_defaults = _make_omc_inputs(atoms)
    calc_defaults |= {"pp_version": "54", "incar_copilot": False}

    return run_and_summarize(
        atoms,
        calc_defaults=calc_defaults,
        calc_swaps=calc_kwargs,
        additional_fields={"name": "OMC Static"} | (additional_fields or {}),
        copy_files=copy_files,
    )


def _make_omc_inputs(atoms: Atoms) -> dict:
    """
    Helper function to make a fairchem input set.

    Parameters
    ----------
    atoms
        Atoms object
    dataset
        Dataset to use. Currently only "omc" is supported.

    Returns
    -------
    dict
        Dictionary of ASE VASP calculator parameters.
    """
    from atomate2.vasp.sets.core import StaticSetGenerator

    input_generator = StaticSetGenerator(
        user_incar_settings={
            "ADDGRID": True,
            "ALGO": "Normal",
            "EDIFF": 1e-06,
            "ENCUT": 520,
            "GGA": "PE",
            "IBRION": -1,
            "ISIF": 0,
            "ISMEAR": 0,
            "ISPIN": 1,
            "IVDW": 11,
            "LREAL": False,
            "LMIXTAU": True,
            "LASPH": True,
            "LORBIT": 11,
            "LWAVE": False,
            "LAECHG": False,
            "LVTOT": False,
            "NELM": 200,
            "NELMDL": -10,
            "NSW": 0,
            "PREC": "Normal",
            "SIGMA": 0.1,
        },
        user_potcar_functional="PBE_54_W_HASH",
        auto_kspacing=True,
    )
    return MPtoASEConverter(atoms=atoms).convert_input_generator(input_generator)

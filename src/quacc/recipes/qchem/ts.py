"""Transition state recipes for Q-Chem."""

from __future__ import annotations

from typing import TYPE_CHECKING

from monty.dev import requires

from quacc import SETTINGS, job, strip_decorator
from quacc.recipes.qchem._base import base_opt_fn
from quacc.recipes.qchem.core import _BASE_SET, relax_job
from quacc.utils.dicts import recursive_dict_merge

try:
    from sella import IRC, Sella

except ImportError:
    Sella = False

if TYPE_CHECKING:
    from typing import Any, Literal

    from ase.atoms import Atoms

    from quacc.schemas._aliases.ase import OptSchema
    from quacc.utils.files import Filenames, SourceDirectory


@job
@requires(Sella, "Sella must be installed. Refer to the quacc documentation.")
def ts_job(
    atoms: Atoms,
    charge: int,
    spin_multiplicity: int,
    method: str = "wb97mv",
    basis: str = "def2-svpd",
    opt_params: dict[str, Any] | None = None,
    copy_files: SourceDirectory | dict[SourceDirectory, Filenames] | None = None,
    **calc_kwargs,
) -> OptSchema:
    """
    TS optimize a molecular structure.

    Parameters
    ----------
    atoms
        Atoms object
    charge
        Charge of the system.
    spin_multiplicity
        Multiplicity of the system.
    method
        DFT exchange-correlation functional or other electronic structure
        method. Defaults to wB97M-V.
    basis
        Basis set. Defaults to def2-SVPD.
    opt_params
        Dictionary of custom kwargs for the optimization process. Set a value
        to `quacc.Remove` to remove a pre-existing key entirely. For a list of available
        keys, refer to [quacc.runners.ase.run_opt][].
    copy_files
        Files to copy (and decompress) from source to the runtime directory.
    **calc_kwargs
        Custom kwargs for the calculator. Set a value to `quacc.Remove` to remove
        a pre-existing key entirely. See [quacc.calculators.qchem.qchem.QChem][] for more
        details.

    Returns
    -------
    OptSchema
        Dictionary of results from [quacc.schemas.ase.summarize_opt_run][].
        See the type-hint for the data structure.
    """

    calc_defaults = recursive_dict_merge(
        _BASE_SET, {"rem": {"job_type": "force", "method": method, "basis": basis}}
    )
    opt_defaults = {"optimizer": Sella, "optimizer_kwargs": {"order": 1}}

    if opt_params and opt_params.get("optimizer", Sella) is not Sella:
        raise ValueError("Only Sella should be used for TS optimization.")

    return base_opt_fn(
        atoms,
        charge,
        spin_multiplicity,
        calc_defaults=calc_defaults,
        calc_swaps=calc_kwargs,
        opt_defaults=opt_defaults,
        opt_params=opt_params,
        additional_fields={"name": "Q-Chem TS"},
        copy_files=copy_files,
    )


@job
@requires(Sella, "Sella must be installed. Refer to the quacc documentation.")
def irc_job(
    atoms: Atoms,
    charge: int,
    spin_multiplicity: int,
    direction: Literal["forward", "reverse"] = "forward",
    method: str = "wb97mv",
    basis: str = "def2-svpd",
    opt_params: dict[str, Any] | None = None,
    copy_files: SourceDirectory | dict[SourceDirectory, Filenames] | None = None,
    **calc_kwargs,
) -> OptSchema:
    """
    IRC optimize a molecular structure.

    Parameters
    ----------
    atoms
        Atoms object
    charge
        Charge of the system.
    spin_multiplicity
        Multiplicity of the system.
    direction
        Direction of the IRC. Should be "forward" or "reverse".
    method
        DFT exchange-correlation functional or other electronic structure
        method. Defaults to wB97M-V.
    basis
        Basis set. Defaults to def2-SVPD.
    opt_params
        Dictionary of custom kwargs for the optimization process. Set a value
        to `quacc.Remove` to remove a pre-existing key entirely. For a list of available
        keys, refer to [quacc.runners.ase.run_opt][].
    copy_files
        Files to copy (and decompress) from source to the runtime directory.
    **calc_kwargs
        Custom kwargs for the calculator. Set a value to `quacc.Remove` to remove
        a pre-existing key entirely. See [quacc.calculators.qchem.qchem.QChem][] for more
        details.

    Returns
    -------
    OptSchema
        Dictionary of results from [quacc.schemas.ase.summarize_opt_run][]
    """

    calc_defaults = recursive_dict_merge(
        _BASE_SET, {"rem": {"job_type": "force", "method": method, "basis": basis}}
    )
    opt_defaults = {
        "optimizer": IRC,
        "optimizer_kwargs": {"keep_going": True},
        "run_kwargs": {"direction": direction},
    }
    if opt_params and opt_params.get("optimizer", IRC) is not IRC:
        raise ValueError("Only Sella's IRC should be used for IRC optimization.")

    return base_opt_fn(
        atoms,
        charge,
        spin_multiplicity,
        calc_defaults=calc_defaults,
        calc_swaps=calc_kwargs,
        opt_defaults=opt_defaults,
        opt_params=opt_params,
        additional_fields={"name": "Q-Chem IRC"},
        copy_files=copy_files,
    )


@job
@requires(Sella, "Sella must be installed. Refer to the quacc documentation.")
def quasi_irc_job(
    atoms: Atoms,
    charge: int,
    spin_multiplicity: int,
    direction: Literal["forward", "reverse"] = "forward",
    method: str = "wb97mv",
    basis: str = "def2-svpd",
    irc_job_kwargs: dict[str, Any] | None = None,
    relax_job_kwargs: dict[str, Any] | None = None,
    copy_files: SourceDirectory | dict[SourceDirectory, Filenames] | None = None,
) -> OptSchema:
    """
    Quasi-IRC optimize a molecular structure. Runs `irc_job` for 10 steps (default)
    followed by `relax_job`.

    Parameters
    ----------
    atoms
        Atoms object.
    charge
        Charge of the system.
    spin_multiplicity
        Multiplicity of the system.
    direction
        Direction of the IRC. Should be "forward" or "reverse".
    irc_job_kwargs
        Dictionary of kwargs for the `irc_job`.
    relax_job_kwargs
        Dictionary of kwargs for the `relax_job`.
    copy_files
        Files to copy (and decompress) from source to the runtime directory.

    Returns
    -------
    OptSchema
        Dictionary of results from [quacc.schemas.ase.summarize_opt_run][]
    """

    default_settings = SETTINGS.model_copy()

    irc_job_defaults = {
        "charge": charge,
        "spin_multiplicity": spin_multiplicity,
        "direction": direction,
        "method": method,
        "basis": basis,
        "opt_params": {"max_steps": 10},
        "copy_files": copy_files,
    }
    relax_job_defaults = {
        "charge": charge,
        "spin_multiplicity": spin_multiplicity,
        "method": method,
        "basis": basis,
    }
    irc_job_kwargs = recursive_dict_merge(irc_job_defaults, irc_job_kwargs)
    relax_job_kwargs = recursive_dict_merge(relax_job_defaults, relax_job_kwargs)

    SETTINGS.CHECK_CONVERGENCE = False
    irc_summary = strip_decorator(irc_job)(atoms, **irc_job_kwargs)

    SETTINGS.CHECK_CONVERGENCE = default_settings.CHECK_CONVERGENCE
    relax_summary = strip_decorator(relax_job)(irc_summary["atoms"], **relax_job_kwargs)

    relax_summary["initial_irc"] = irc_summary

    return relax_summary

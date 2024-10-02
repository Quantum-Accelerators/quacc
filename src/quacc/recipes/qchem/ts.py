"""Transition state recipes for Q-Chem."""

from __future__ import annotations

from importlib.util import find_spec
from typing import TYPE_CHECKING

from monty.dev import requires

from quacc import job
from quacc.atoms.core import perturb
from quacc.recipes.qchem._base import run_and_summarize_opt
from quacc.recipes.qchem.core import _BASE_SET
from quacc.utils.dicts import recursive_dict_merge

has_sella = bool(find_spec("sella"))
if has_sella:
    from sella import IRC, Sella

if TYPE_CHECKING:
    from typing import Any, Literal

    from ase.atoms import Atoms
    from numpy.typing import NDArray

    from quacc.types import Filenames, OptParams, OptSchema, SourceDirectory


@job
@requires(has_sella, "Sella must be installed. Refer to the quacc documentation.")
def ts_job(
    atoms: Atoms,
    charge: int = 0,
    spin_multiplicity: int = 1,
    method: str = "wb97mv",
    basis: str = "def2-svpd",
    opt_params: OptParams | None = None,
    copy_files: SourceDirectory | dict[SourceDirectory, Filenames] | None = None,
    additional_fields: dict[str, Any] | None = None,
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
        method.
    basis
        Basis set.
    opt_params
        Dictionary of custom kwargs for the optimization process. For a list
        of available keys, refer to [quacc.runners.ase.Runner.run_opt][].
    copy_files
        Files to copy (and decompress) from source to the runtime directory.
    **calc_kwargs
        Custom kwargs for the calculator. Set a value to `quacc.Remove` to remove
        a pre-existing key entirely. See [quacc.calculators.qchem.qchem.QChem][] for more
        details.

    Returns
    -------
    OptSchema
        Dictionary of results from [quacc.schemas.ase.Summarize.opt][].
        See the type-hint for the data structure.
    """
    calc_defaults = recursive_dict_merge(
        _BASE_SET, {"rem": {"job_type": "force", "method": method, "basis": basis}}
    )
    opt_defaults = {"optimizer": Sella, "optimizer_kwargs": {"order": 1}}

    if opt_params and opt_params.get("optimizer", Sella) is not Sella:
        raise ValueError("Only Sella should be used for TS optimization.")

    return run_and_summarize_opt(
        atoms,
        charge,
        spin_multiplicity,
        calc_defaults=calc_defaults,
        calc_swaps=calc_kwargs,
        opt_defaults=opt_defaults,
        opt_params=opt_params,
        additional_fields={"name": "Q-Chem TS"} | (additional_fields or {}),
        copy_files=copy_files,
    )


@job
@requires(has_sella, "Sella must be installed. Refer to the quacc documentation.")
def irc_job(
    atoms: Atoms,
    charge: int = 0,
    spin_multiplicity: int = 1,
    direction: Literal["forward", "reverse"] = "forward",
    method: str = "wb97mv",
    basis: str = "def2-svpd",
    opt_params: OptParams | None = None,
    copy_files: SourceDirectory | dict[SourceDirectory, Filenames] | None = None,
    additional_fields: dict[str, Any] | None = None,
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
        method.
    basis
        Basis set.
    opt_params
        Dictionary of custom kwargs for the optimization process. For a list
        of available keys, refer to [quacc.runners.ase.Runner.run_opt][].
    copy_files
        Files to copy (and decompress) from source to the runtime directory.
    additional_fields
        Additional fields to add to the results dictionary.
    **calc_kwargs
        Custom kwargs for the calculator. Set a value to `quacc.Remove` to remove
        a pre-existing key entirely. See [quacc.calculators.qchem.qchem.QChem][] for more
        details.

    Returns
    -------
    OptSchema
        Dictionary of results from [quacc.schemas.ase.Summarize.opt][]
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

    return run_and_summarize_opt(
        atoms,
        charge,
        spin_multiplicity,
        calc_defaults=calc_defaults,
        calc_swaps=calc_kwargs,
        opt_defaults=opt_defaults,
        opt_params=opt_params,
        additional_fields={"name": "Q-Chem IRC"} | (additional_fields or {}),
        copy_files=copy_files,
    )


@job
@requires(has_sella, "Sella must be installed. Refer to the quacc documentation.")
def quasi_irc_job(
    atoms: Atoms,
    mode: list[list[float]] | NDArray,
    perturb_magnitude: float = 0.6,
    direction: Literal["forward", "reverse"] = "forward",
    charge: int = 0,
    spin_multiplicity: int = 1,
    method: str = "wb97mv",
    basis: str = "def2-svpd",
    opt_params: OptParams | None = None,
    copy_files: SourceDirectory | dict[SourceDirectory, Filenames] | None = None,
    additional_fields: dict[str, Any] | None = None,
    **calc_kwargs,
) -> OptSchema:
    """
    Quasi-IRC to optimize a reaction endpoint from a transition-state with known vibrational frequency modes.
    Perturbs the structure of `atoms` by a finite amount (0.6 * the normalized mode magnitude) along the specified
    vibrational frequency mode (assumed to be the transition mode), and then performs a `relax_job` on the perturbed
    structure.

    Parameters
    ----------
    atoms
        Atoms object.
    mode
        Transition mode. This should be an Nx3 matrix, where N is the number of atoms in `atoms`.
    perturb_magnitude
        Factor to multiply the transition mode. Default is 0.6. In some cases, it may be advisable to increase this
        factor, perhaps to 1.0 or 1.1. Lowering it is not generally found to be helpful.
    direction
        Direction of the (Quasi)IRC. Should be "forward" or "reverse".
    charge
        Charge of the system.
    spin_multiplicity
        Multiplicity of the system.
    method
        DFT exchange-correlation functional or other electronic structure
        method.
    basis
        Basis set.
    opt_params
        Dictionary of custom kwargs for the optimization process. For a list
        of available keys, refer to [quacc.runners.ase.Runner.run_opt][].
    copy_files
        Files to copy (and decompress) from source to the runtime directory.
    additional_fields
        Additional fields to add to the results dictionary.
    **calc_kwargs
        Custom kwargs for the calculator. Set a value to `quacc.Remove` to remove
        a pre-existing key entirely. See [quacc.calculators.qchem.qchem.QChem][] for more
        details.


    Returns
    -------
    OptSchema
        Dictionary of results from [quacc.schemas.ase.Summarize.opt][]
    """
    calc_defaults = recursive_dict_merge(
        _BASE_SET, {"rem": {"job_type": "force", "method": method, "basis": basis}}
    )
    opt_defaults = {"optimizer": Sella} if has_sella else {}

    scale = perturb_magnitude if direction == "forward" else perturb_magnitude * -1

    return run_and_summarize_opt(
        perturb(atoms, mode, scale),
        charge=charge,
        spin_multiplicity=spin_multiplicity,
        calc_defaults=calc_defaults,
        calc_swaps=calc_kwargs,
        opt_defaults=opt_defaults,
        opt_params=opt_params,
        additional_fields={"name": "Q-Chem Quasi-IRC"} | (additional_fields or {}),
        copy_files=copy_files,
    )

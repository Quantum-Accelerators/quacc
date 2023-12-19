"""Core recipes for the tblite code."""
from __future__ import annotations

from typing import TYPE_CHECKING

from ase.optimize import FIRE
from monty.dev import requires

from quacc import job
from quacc.runners.ase import run_calc, run_opt
from quacc.schemas.ase import summarize_opt_run, summarize_run, summarize_vib_and_thermo
from quacc.utils.dicts import merge_dicts

try:
    from shutil import which

    from xtb_ase.calculator import XTB, XTBProfile

    has_xtb = which("xtb")
except ImportError:
    has_xtb = False

if TYPE_CHECKING:
    from typing import Any, Literal

    from ase.atoms import Atoms

    from quacc.runners.ase import VibKwargs
    from quacc.schemas._aliases.ase import OptSchema, RunSchema, VibThermoSchema


@job
@requires(
    has_xtb, "xtb and xtb-ase must be installed. Refer to the quacc documentation."
)
def static_job(
    atoms: Atoms,
    method: Literal["gfn0-xtb", "gfn1-xtb", "gfn2-xTB", "gfn-ff"] = "gfn2-xtb",
    charge: int = 0,
    spin_multiplicity: int = 1,
    spinpol: bool | None = None,
    **calc_kwargs,
) -> RunSchema:
    """
    Carry out a single-point calculation.

    Parameters
    ----------
    atoms
        Atoms object
    method
        The xTB method to use. Case-insensitive.
    charge
        The net charge of the system.
    spin_multiplicity
        The spin multiplicity of the system.
    spinpol
        Whether to use spin-polarized xTB. If None, `spinpol` will be automatically
        set to True if `spin_multiplicity` > 1.
    **calc_kwargs
        Any additional xTB parameters to be written out to a detailed input file, e.g. in the format
        of `scc={"temp": 500}`. See https://github.com/grimme-lab/xtb/blob/main/man/xcontrol.7.adoc.

    Returns
    -------
    RunSchema
        Dictionary of results from [quacc.schemas.ase.summarize_run][]
    """

    atoms.calc = XTB(
        method=method,
        charge=charge,
        uhf=spin_multiplicity - 1,
        spinpol=spinpol,
        **calc_kwargs,
    )

    final_atoms = run_calc(atoms)
    return summarize_run(
        final_atoms, input_atoms=atoms, additional_fields={"name": "xTB Static"}
    )


@job
@requires(
    has_xtb, "xtb and xtb-ase must be installed. Refer to the quacc documentation."
)
def relax_job(
    atoms: Atoms,
    method: Literal["gfn0-xtb", "gfn1-xtb", "gfn2-xTB", "gfn-ff"] = "gfn2-xtb",
    charge: int = 0,
    spin_multiplicity: int = 1,
    spinpol: bool | None = None,
    **calc_kwargs,
) -> RunSchema:
    """
    Carry out a relaxation using xTB optimizers.

    Parameters
    ----------
    atoms
        Atoms object
    method
        The xTB method to use. Case-insensitive.
    charge
        The net charge of the system.
    spin_multiplicity
        The spin multiplicity of the system.
    spinpol
        Whether to use spin-polarized xTB. If None, `spinpol` will be automatically
        set to True if `spin_multiplicity` > 1.
    **calc_kwargs
        Any additional xTB parameters to be written out to a detailed input file, e.g. in the format
        of `scc={"temp": 500}`. See https://github.com/grimme-lab/xtb/blob/main/man/xcontrol.7.adoc.

    Returns
    -------
    RunSchema
        Dictionary of results from [quacc.schemas.ase.summarize_run][]
    """
    atoms.calc = XTB(
        profile=XTBProfile(["xtb", "--opt"]),
        method=method,
        charge=charge,
        uhf=spin_multiplicity - 1,
        spinpol=spinpol,
        **calc_kwargs,
    )

    final_atoms = run_calc(atoms)
    return summarize_run(
        final_atoms, input_atoms=atoms, additional_fields={"name": "xTB Static"}
    )


@job
@requires(has_xtb, "xtb-ase must be installed. Refer to the quacc documentation.")
def ase_relax_job(
    atoms: Atoms,
    method: Literal["gfn0-xtb", "gfn1-xtb", "gfn2-xTB", "gfn-ff"] = "gfn2-xtb",
    charge: int = 0,
    spin_multiplicity: int = 1,
    spinpol: bool | None = None,
    opt_params: dict[str, Any] | None = None,
    **calc_kwargs,
) -> OptSchema:
    """
    Carry out a relaxation using ASE optimizers.

    Parameters
    ----------
    atoms
        Atoms object
    method
        The xTB method to use. Case-insensitive.
    charge
        The net charge of the system.
    spin_multiplicity
        The spin multiplicity of the system.
    spinpol
        Whether to use spin-polarized xTB. If None, `spinpol` will be automatically
        set to True if `spin_multiplicity` > 1.
    opt_params
        Dictionary of custom kwargs for the optimization process. Set a value
        to `None` to remove a pre-existing key entirely. For a list of available
        keys, refer to [quacc.runners.ase.run_opt][].
    **calc_kwargs
        Any additional xTB parameters to be written out to a detailed input file, e.g. in the format
        of `scc={"temp": 500}`. See https://github.com/grimme-lab/xtb/blob/main/man/xcontrol.7.adoc.

    Returns
    -------
    OptSchema
        Dictionary of results from [quacc.schemas.ase.summarize_opt_run][]
    """

    atoms.calc = XTB(
        method=method,
        charge=charge,
        uhf=spin_multiplicity - 1,
        spinpol=spinpol,
        **calc_kwargs,
    )
    if atoms.pbc.any():
        raise ValueError("Periodic boundary conditions are not currently supported.")

    opt_defaults = {"fmax": 0.01, "max_steps": 1000, "optimizer": FIRE}
    opt_flags = merge_dicts(opt_defaults, opt_params)

    dyn = run_opt(atoms, **opt_flags)

    return summarize_opt_run(dyn, additional_fields={"name": "xTB Relax"})

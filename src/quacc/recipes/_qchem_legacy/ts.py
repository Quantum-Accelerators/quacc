"""Transition state recipes for the Q-Chem."""
from __future__ import annotations

from typing import TYPE_CHECKING

import psutil
from monty.dev import requires

from quacc import SETTINGS, job, strip_decorator
from quacc.recipes._qchem_legacy._base import base_opt_fn
from quacc.recipes._qchem_legacy.core import relax_job
from quacc.utils.dicts import recursive_dict_merge

try:
    from sella import IRC, Sella

except ImportError:
    Sella = False

if TYPE_CHECKING:
    from pathlib import Path
    from typing import Any, Literal

    from ase.atoms import Atoms

    from quacc.schemas._aliases.ase import OptSchema


@job
@requires(Sella, "Sella must be installed. Refer to the quacc documentation.")
def ts_job(
    atoms: Atoms,
    charge: int,
    spin_multiplicity: int,
    method: str = "wb97mv",
    basis: str = "def2-svpd",
    scf_algorithm: str = "diis",
    pcm_dielectric: str | None = None,
    smd_solvent: str | None = None,
    n_cores: int | None = None,
    overwrite_inputs: dict[str, Any] | None = None,
    opt_params: dict[str, Any] | None = None,
    copy_files: str | Path | list[str | Path] | None = None,
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
    scf_algorithm
        Algorithm used to converge the SCF. Defaults to "diis", but for
        particularly difficult cases, "gdm" should be employed instead.
    pcm_dielectric
        Dielectric constant of the optional polarizable continuum implicit
        solvation model. Defaults to None, in which case PCM will not be
        employed.
    smd_solvent
        Solvent to use for SMD implicit solvation model. Examples include
        "water", "ethanol", "methanol", and "acetonitrile". Refer to the Q-Chem
        manual for a complete list of solvents available. Defaults to None, in
        which case SMD will not be employed.
    n_cores
        Number of cores to use for the Q-Chem calculation. Defaults to use all
        cores available on a given node.
    overwrite_inputs
        Dictionary passed to `pymatgen.io.qchem.QChemDictSet` which can modify
        default values set therein as well as set additional Q-Chem parameters.
        See QChemDictSet documentation for more details.
    opt_params
        Dictionary of custom kwargs for [quacc.runners.ase.run_opt][]
    copy_files
        File(s) to copy to the runtime directory. If a directory is provided, it will be recursively unpacked.

    Returns
    -------
    OptSchema
        Dictionary of results from [quacc.schemas.ase.summarize_opt_run][]
    """

    calc_defaults = {
        "basis_set": basis,
        "scf_algorithm": scf_algorithm,
        "method": method,
        "charge": charge,
        "spin_multiplicity": spin_multiplicity,
        "cores": n_cores or psutil.cpu_count(logical=False),
        "qchem_input_params": {
            "pcm_dielectric": pcm_dielectric,
            "smd_solvent": smd_solvent,
            "overwrite_inputs": overwrite_inputs,
            "max_scf_cycles": 200 if scf_algorithm.lower() == "gdm" else 100,
        },
    }
    opt_defaults = {
        "fmax": 0.01,
        "max_steps": 1000,
        "optimizer": Sella,
        "optimizer_kwargs": {"order": 1},
    }

    if opt_params and opt_params.get("optimizer", Sella) is not Sella:
        raise ValueError("Only Sella should be used for TS optimization.")

    return base_opt_fn(
        atoms,
        charge,
        spin_multiplicity,
        calc_defaults=calc_defaults,
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
    scf_algorithm: str = "diis",
    pcm_dielectric: str | None = None,
    smd_solvent: str | None = None,
    n_cores: int | None = None,
    overwrite_inputs: dict[str, Any] | None = None,
    opt_params: dict[str, Any] | None = None,
    copy_files: str | Path | list[str | Path] | None = None,
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
    scf_algorithm
        Algorithm used to converge the SCF. Defaults to "diis", but for
        particularly difficult cases, "gdm" should be employed instead.
    pcm_dielectric
        Dielectric constant of the optional polarizable continuum implicit
        solvation model. Defaults to None, in which case PCM will not be
        employed.
    smd_solvent
        Solvent to use for SMD implicit solvation model. Examples include
        "water", "ethanol", "methanol", and "acetonitrile". Refer to the Q-Chem
        manual for a complete list of solvents available. Defaults to None, in
        which case SMD will not be employed.
    n_cores
        Number of cores to use for the Q-Chem calculation. Defaults to use all
        cores available on a given node.
    overwrite_inputs
        Dictionary passed to `pymatgen.io.qchem.QChemDictSet` which can modify
        default values set therein as well as set additional Q-Chem parameters.
        See QChemDictSet documentation for more details.
    opt_params
        Dictionary of custom kwargs for [quacc.runners.ase.run_opt][]
    copy_files
        File(s) to copy to the runtime directory. If a directory is provided, it will be recursively unpacked.

    Returns
    -------
    OptSchema
        Dictionary of results from [quacc.schemas.ase.summarize_opt_run][]
    """

    calc_defaults = {
        "basis_set": basis,
        "scf_algorithm": scf_algorithm,
        "method": method,
        "charge": charge,
        "spin_multiplicity": spin_multiplicity,
        "cores": n_cores or psutil.cpu_count(logical=False),
        "qchem_input_params": {
            "pcm_dielectric": pcm_dielectric,
            "smd_solvent": smd_solvent,
            "overwrite_inputs": overwrite_inputs,
            "max_scf_cycles": 200 if scf_algorithm.lower() == "gdm" else 100,
        },
    }
    opt_defaults = {
        "fmax": 0.01,
        "max_steps": 1000,
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
    scf_algorithm: str = "diis",
    pcm_dielectric: str | None = None,
    smd_solvent: str | None = None,
    n_cores: int | None = None,
    overwrite_inputs: dict[str, Any] | None = None,
    irc_opt_params: dict[str, Any] | None = None,
    relax_opt_params: dict[str, Any] | None = None,
    copy_files: str | Path | list[str | Path] | None = None,
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
    irc_opt_params
        Dictionary of opt_params kwargs for the irc_job.
    relax_opt_params
        Dictionary of opt_params kwargs for the relax_job.
    copy_files
        File(s) to copy to the runtime directory. If a directory is provided, it will be recursively unpacked.

    Returns
    -------
    OptSchema
        Dictionary of results from [quacc.schemas.ase.summarize_opt_run][]
    """

    default_settings = SETTINGS.model_copy()

    irc_opt_params_defaults = {"fmax": 100, "max_steps": 10}
    irc_opt_params = recursive_dict_merge(irc_opt_params_defaults, irc_opt_params)

    SETTINGS.CHECK_CONVERGENCE = False
    irc_summary = strip_decorator(irc_job)(
        atoms,
        charge,
        spin_multiplicity,
        direction=direction,
        method=method,
        basis=basis,
        scf_algorithm=scf_algorithm,
        pcm_dielectric=pcm_dielectric,
        smd_solvent=smd_solvent,
        n_cores=n_cores,
        overwrite_inputs=overwrite_inputs,
        opt_params=irc_opt_params,
        copy_files=copy_files,
    )

    SETTINGS.CHECK_CONVERGENCE = default_settings.CHECK_CONVERGENCE
    relax_summary = strip_decorator(relax_job)(
        irc_summary["atoms"],
        charge,
        spin_multiplicity,
        method=method,
        basis=basis,
        scf_algorithm=scf_algorithm,
        pcm_dielectric=pcm_dielectric,
        smd_solvent=smd_solvent,
        n_cores=n_cores,
        overwrite_inputs=overwrite_inputs,
        opt_params=relax_opt_params,
    )

    relax_summary["initial_irc"] = irc_summary

    return relax_summary

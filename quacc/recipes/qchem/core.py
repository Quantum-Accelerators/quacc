"""
Core recipes for the Q-Chem
"""
from __future__ import annotations

from copy import deepcopy
from typing import Literal

import covalent as ct
from ase.atoms import Atoms
from monty.dev import requires

from quacc.calculators.qchem import QChem
from quacc.schemas.ase import summarize_opt_run, summarize_run
from quacc.util.calc import run_ase_opt, run_calc

try:
    from sella import Sella, IRC
except ImportError:
    Sella = None


@ct.electron
def static_job(
    atoms: Atoms,
    cores: int | None = None,
    charge: int | None = None,
    mult: int | None = None,
    xc: str = "wb97mv",
    basis: str = "def2-tzvpd",
    scf_algorithm: str = "diis",
    pcm_dielectric: str | None = None,
    smd_solvent: str | None = None,
    swaps: dict | None = None,
) -> dict:
    """
    Carry out a single-point calculation.

    Parameters
    ----------
    atoms
        Atoms object.
    cores
        Number of cores to use for the Q-Chem calculation.
        Effectively defaults to use all cores available on a given node, so this only needs to
        be set by the user if less than all available cores should be used.
    charge
        The total charge of the molecular system.
        Effectively defaults to zero.
    mult
        The spin multiplicity of the molecular system.
        Effectively defaults to the lowest spin state given the molecular structure and charge.
    xc
        Exchange-correlation functional.
        Defaults to wB97M-V.
    basis
        Basis set.
        Defaults to def2-TZVPD.
    scf_algorithm
        Algorithm used to converge the SCF.
        Defaults to DIIS, but for particularly difficult cases, GDM should be employed instead.
    pcm_dielectric
        Dielectric constant of the optional polarizable continuum impicit solvation model.
        Defaults no None, in which case PCM will not be employed.
    smd_solvent
        Solvent to use for SMD implicit solvation model. Examples include "water", "ethanol", "methanol",
        and "acetonitrile". Refer to the Q-Chem manual for a complete list of solvents available.
        Defulats to None, in which case SMD will not be employed.
    swaps
        Dictionary of custom kwargs for the calculator. Must be formatted consistently with Pymatgen's
        QChemDictSet's overwrite_inputs. For example: {"rem": {"symmetry": "true"}}.

    Returns
    -------
    dict
        Dictionary of results from quacc.schemas.ase.summarize_run
    """

    if pcm_dielectric is not None and smd_solvent is not None:
        raise RuntimeError("PCM and SMD cannot be employed simultaneously! Exiting...")

    input_atoms = deepcopy(atoms)

    overwrite_inputs = {"rem": {"method": xc}}

    swaps = swaps or {}
    for key in swaps:
        if key not in overwrite_inputs:
            overwrite_inputs[key] = {}
            for subkey in swaps[key]:
                overwrite_inputs[key][subkey] = swaps[key][subkey]

    qchem_input_params = {
        "basis_set": basis,
        "scf_algorithm": scf_algorithm,
        "qchem_version": 6,
        "pcm_dielectric": pcm_dielectric,
        "smd_solvent": smd_solvent,
        "overwrite_inputs": overwrite_inputs,
    }

    if scf_algorithm.lower() == "gdm":
        qchem_input_params["max_scf_cycles"] = 200

    calc = QChem(
        input_atoms=atoms,
        cores=cores,
        charge=charge,
        spin_multiplicity=mult,
        qchem_input_params=qchem_input_params,
    )
    atoms.calc = calc
    atoms = run_calc(atoms)
    return summarize_run(
        atoms,
        input_atoms=input_atoms,
        additional_fields={"name": "Q-Chem Static"},
    )


@ct.electron
@requires(
    Sella,
    "Sella must be installed. pip install sella",
)
def relax_job(
    atoms: Atoms,
    cores: int | None = None,
    charge: int | None = None,
    mult: int | None = None,
    xc: str = "wb97mv",
    basis: str = "def2-svpd",
    scf_algorithm: str = "diis",
    pcm_dielectric: str | None = None,
    smd_solvent: str | None = None,
    swaps: dict | None = None,
    opt_swaps: dict | None = None,
    check_convergence: bool = True,
) -> dict:
    """
    Optimize aka "relax" a molecular structure.

    Parameters
    ----------
    atoms
        Atoms object.
    cores
        Number of cores to use for the Q-Chem calculation.
        Effectively defaults to use all cores available on a given node, so this only needs to
        be set by the user if less than all available cores should be used.
    charge
        The total charge of the molecular system.
        Effectively defaults to zero.
    mult
        The spin multiplicity of the molecular system.
        Effectively defaults to the lowest spin state given the molecular structure and charge.
    xc
        Exchange-correlation functional.
        Defaults to wB97M-V.
    basis
        Basis set.
        Defaults to def2-SVPD.
    scf_algorithm
        Algorithm used to converge the SCF.
        Defaults to DIIS, but for particularly difficult cases, GDM should be employed instead.
    pcm_dielectric
        Dielectric constant of the optional polarizable continuum impicit solvation model.
        Defaults no None, in which case PCM will not be employed.
    smd_solvent
        Solvent to use for SMD implicit solvation model. Examples include "water", "ethanol", "methanol",
        and "acetonitrile". Refer to the Q-Chem manual for a complete list of solvents available.
        Defulats to None, in which case SMD will not be employed.
    swaps
        Dictionary of custom kwargs for the calculator. Must be formatted consistently with Pymatgen's
        QChemDictSet's overwrite_inputs. For example: {"rem": {"symmetry": "true"}}.
    opt_swaps
        Dictionary of custom kwargs for run_ase_opt
            opt_defaults = {"fmax": 0.01, "max_steps": 1000, "optimizer": "Sella"}

    Returns
    -------
    dict
        Dictionary of results from quacc.schemas.ase.summarize_opt_run
    """

    # Reminder to self: exposing TRICs?

    opt_swaps = opt_swaps or {}
    opt_defaults = {
        "fmax": 0.01,
        "max_steps": 1000,
        "optimizer": Sella,
        "optimizer_kwargs": {},
    }
    opt_flags = opt_defaults | opt_swaps
    if "sella.optimize" in opt_flags["optimizer"].__module__:
        opt_flags["optimizer_kwargs"]["order"] = 0

    if pcm_dielectric is not None and smd_solvent is not None:
        raise ValueError("PCM and SMD cannot be employed simultaneously! Exiting...")

    overwrite_inputs = {"rem": {"method": xc}}

    swaps = swaps or {}
    for key in swaps:
        if key not in overwrite_inputs:
            overwrite_inputs[key] = {}
            for subkey in swaps[key]:
                overwrite_inputs[key][subkey] = swaps[key][subkey]

    qchem_input_params = {
        "basis_set": basis,
        "scf_algorithm": scf_algorithm,
        "qchem_version": 6,
        "pcm_dielectric": pcm_dielectric,
        "smd_solvent": smd_solvent,
        "overwrite_inputs": overwrite_inputs,
    }

    if scf_algorithm.lower() == "gdm":
        qchem_input_params["max_scf_cycles"] = 200

    calc = QChem(
        input_atoms=atoms,
        cores=cores,
        charge=charge,
        spin_multiplicity=mult,
        qchem_input_params=qchem_input_params,
    )
    atoms.calc = calc
    dyn = run_ase_opt(atoms, **opt_flags)

    return summarize_opt_run(
        dyn,
        check_convergence=check_convergence,
        additional_fields={"name": "Q-Chem Optimization"},
    )


@ct.electron
@requires(
    Sella,
    "Sella must be installed. pip install sella",
)
def ts_job(
    atoms: Atoms,
    cores: int | None = None,
    charge: int | None = None,
    mult: int | None = None,
    xc: str = "wb97mv",
    basis: str = "def2-svpd",
    scf_algorithm: str = "diis",
    pcm_dielectric: str | None = None,
    smd_solvent: str | None = None,
    swaps: dict | None = None,
    opt_swaps: dict | None = None,
    check_convergence: bool = True,
) -> dict:
    """
    TS optimize a molecular structure.

    Parameters
    ----------
    atoms
        Atoms object.
    cores
        Number of cores to use for the Q-Chem calculation.
        Effectively defaults to use all cores available on a given node, so this only needs to
        be set by the user if less than all available cores should be used.
    charge
        The total charge of the molecular system.
        Effectively defaults to zero.
    mult
        The spin multiplicity of the molecular system.
        Effectively defaults to the lowest spin state given the molecular structure and charge.
    xc
        Exchange-correlation functional.
        Defaults to wB97M-V.
    basis
        Basis set.
        Defaults to def2-SVPD.
    scf_algorithm
        Algorithm used to converge the SCF.
        Defaults to DIIS, but for particularly difficult cases, GDM should be employed instead.
    pcm_dielectric
        Dielectric constant of the optional polarizable continuum impicit solvation model.
        Defaults no None, in which case PCM will not be employed.
    smd_solvent
        Solvent to use for SMD implicit solvation model. Examples include "water", "ethanol", "methanol",
        and "acetonitrile". Refer to the Q-Chem manual for a complete list of solvents available.
        Defulats to None, in which case SMD will not be employed.
    swaps
        Dictionary of custom kwargs for the calculator. Must be formatted consistently with Pymatgen's
        QChemDictSet's overwrite_inputs. For example: {"rem": {"symmetry": "true"}}.
    opt_swaps
        Dictionary of custom kwargs for run_ase_opt
            opt_defaults = {"fmax": 0.01, "max_steps": 1000, "optimizer": "Sella"}

    Returns
    -------
    dict
        Dictionary of results from quacc.schemas.ase.summarize_opt_run
    """

    # Reminder to self: exposing TRICs?

    opt_swaps = opt_swaps or {}
    opt_defaults = {
        "fmax": 0.01,
        "max_steps": 1000,
        "optimizer": Sella,
        "optimizer_kwargs": {},
    }
    opt_flags = opt_defaults | opt_swaps
    if "sella.optimize" not in opt_flags["optimizer"].__module__:
        raise RuntimeError("Only Sella should be used for TS optimization! Exiting...")

    if pcm_dielectric is not None and smd_solvent is not None:
        raise ValueError("PCM and SMD cannot be employed simultaneously! Exiting...")

    overwrite_inputs = {"rem": {"method": xc}}

    swaps = swaps or {}
    for key in swaps:
        if key not in overwrite_inputs:
            overwrite_inputs[key] = {}
            for subkey in swaps[key]:
                overwrite_inputs[key][subkey] = swaps[key][subkey]

    qchem_input_params = {
        "basis_set": basis,
        "scf_algorithm": scf_algorithm,
        "qchem_version": 6,
        "pcm_dielectric": pcm_dielectric,
        "smd_solvent": smd_solvent,
        "overwrite_inputs": overwrite_inputs,
    }

    if scf_algorithm.lower() == "gdm":
        qchem_input_params["max_scf_cycles"] = 200

    calc = QChem(
        input_atoms=atoms,
        cores=cores,
        charge=charge,
        spin_multiplicity=mult,
        qchem_input_params=qchem_input_params,
    )
    atoms.calc = calc
    dyn = run_ase_opt(atoms, **opt_flags)

    return summarize_opt_run(
        dyn,
        check_convergence=check_convergence,
        additional_fields={"name": "Q-Chem Optimization"},
    )


@ct.electron
@requires(
    Sella,
    "Sella must be installed. pip install sella",
)
def irc_job(
    atoms: Atoms,
    direction: str,
    cores: int | None = None,
    charge: int | None = None,
    mult: int | None = None,
    xc: str = "wb97mv",
    basis: str = "def2-svpd",
    scf_algorithm: str = "diis",
    pcm_dielectric: str | None = None,
    smd_solvent: str | None = None,
    swaps: dict | None = None,
    opt_swaps: dict | None = None,
    check_convergence: bool = True,
) -> dict:
    """
    TS optimize a molecular structure.

    Parameters
    ----------
    atoms
        Atoms object.
    direction
        Direction of the IRC. Should be "forward" or "reverse".
    cores
        Number of cores to use for the Q-Chem calculation.
        Effectively defaults to use all cores available on a given node, so this only needs to
        be set by the user if less than all available cores should be used.
    charge
        The total charge of the molecular system.
        Effectively defaults to zero.
    mult
        The spin multiplicity of the molecular system.
        Effectively defaults to the lowest spin state given the molecular structure and charge.
    xc
        Exchange-correlation functional.
        Defaults to wB97M-V.
    basis
        Basis set.
        Defaults to def2-SVPD.
    scf_algorithm
        Algorithm used to converge the SCF.
        Defaults to DIIS, but for particularly difficult cases, GDM should be employed instead.
    pcm_dielectric
        Dielectric constant of the optional polarizable continuum impicit solvation model.
        Defaults no None, in which case PCM will not be employed.
    smd_solvent
        Solvent to use for SMD implicit solvation model. Examples include "water", "ethanol", "methanol",
        and "acetonitrile". Refer to the Q-Chem manual for a complete list of solvents available.
        Defulats to None, in which case SMD will not be employed.
    swaps
        Dictionary of custom kwargs for the calculator. Must be formatted consistently with Pymatgen's
        QChemDictSet's overwrite_inputs. For example: {"rem": {"symmetry": "true"}}.
    opt_swaps
        Dictionary of custom kwargs for run_ase_opt
            opt_defaults = {"fmax": 0.01, "max_steps": 1000, "optimizer": "Sella"}

    Returns
    -------
    dict
        Dictionary of results from quacc.schemas.ase.summarize_opt_run
    """

    # Reminder to self: exposing TRICs?

    if direction not in ["forward", "reverse"]:
        raise ValueError("direction must be 'forward' or 'reverse'! Exiting...")

    opt_swaps = opt_swaps or {}
    opt_defaults = {
        "fmax": 0.01,
        "max_steps": 1000,
        "optimizer": IRC,
        "optimizer_kwargs": {},
        "run_kwargs": {"direction": direction}
    }
    opt_flags = opt_defaults | opt_swaps
    if "sella.optimize" not in opt_flags["optimizer"].__module__:
        raise RuntimeError("Only Sella should be used for TS optimization! Exiting...")

    if pcm_dielectric is not None and smd_solvent is not None:
        raise ValueError("PCM and SMD cannot be employed simultaneously! Exiting...")

    overwrite_inputs = {"rem": {"method": xc}}

    swaps = swaps or {}
    for key in swaps:
        if key not in overwrite_inputs:
            overwrite_inputs[key] = {}
            for subkey in swaps[key]:
                overwrite_inputs[key][subkey] = swaps[key][subkey]

    qchem_input_params = {
        "basis_set": basis,
        "scf_algorithm": scf_algorithm,
        "qchem_version": 6,
        "pcm_dielectric": pcm_dielectric,
        "smd_solvent": smd_solvent,
        "overwrite_inputs": overwrite_inputs,
    }

    if scf_algorithm.lower() == "gdm":
        qchem_input_params["max_scf_cycles"] = 200

    calc = QChem(
        input_atoms=atoms,
        cores=cores,
        charge=charge,
        spin_multiplicity=mult,
        qchem_input_params=qchem_input_params,
    )
    atoms.calc = calc
    dyn = run_ase_opt(atoms, **opt_flags)

    return summarize_opt_run(
        dyn,
        check_convergence=check_convergence,
        additional_fields={"name": "Q-Chem Optimization"},
    )

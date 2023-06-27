"""
Core recipes for the Q-Chem
"""
from __future__ import annotations

from copy import deepcopy

import covalent as ct
from ase.atoms import Atoms
from ase.optimize import FIRE
from monty.dev import requires

from quacc.calculators.qchem import QChem
from quacc.schemas.ase import summarize_opt_run, summarize_run
from quacc.util.calc import run_ase_opt, run_calc
from quacc.util.atoms import check_charge_and_spin
from typing import Literal

try:
    from sella import IRC, Sella
except ImportError:
    Sella = None


@ct.electron
def static_job(
    atoms: Atoms,
    cores: int | None = None,
    charge: int | None = None,
    spin_multiplicity: int | None = None,
    method: str = "wb97mv",
    basis: str = "def2-tzvpd",
    scf_algorithm: str = "diis",
    pcm_dielectric: str | None = None,
    smd_solvent: str | None = None,
    overwrite_inputs: dict | None = None,
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
    spin_multiplicity
        The spin multiplicity of the molecular system.
        Effectively defaults to the lowest spin state given the molecular structure and charge.
    method
        DFT exchange-correlation functional or other electronic structure method.
        Defaults to wB97M-V.
    basis
        Basis set.
        Defaults to def2-TZVPD.
    scf_algorithm
        Algorithm used to converge the SCF.
        Defaults to "diis", but for particularly difficult cases, "gdm" should be employed instead.
    pcm_dielectric
        Dielectric constant of the optional polarizable continuum impicit solvation model.
        Defaults to None, in which case PCM will not be employed.
    smd_solvent
        Solvent to use for SMD implicit solvation model. Examples include "water", "ethanol", "methanol",
        and "acetonitrile". Refer to the Q-Chem manual for a complete list of solvents available.
        Defaults to None, in which case SMD will not be employed.
    overwrite_inputs
        Dictionary passed to pymatgen.io.qchem.QChemDictSet which can modify default values set therein
        as well as set additional Q-Chem parameters. See QChemDictSet documentation for more details.

    Returns
    -------
    dict
        Dictionary of results from quacc.schemas.ase.summarize_run
    """

    if pcm_dielectric is not None and smd_solvent is not None:
        raise ValueError("PCM and SMD cannot be employed simultaneously! Exiting...")

    checked_charge, checked_spin_multiplicity = check_charge_and_spin(atoms, charge, spin_multiplicity)

    input_atoms = deepcopy(atoms)

    overwrite_inputs = overwrite_inputs or {}
    if "rem" not in overwrite_inputs:
        overwrite_inputs["rem"] = {}
    if "method" not in overwrite_inputs["rem"]:
        overwrite_inputs["rem"]["method"] = method

    qchem_input_params = {
        "basis_set": basis,
        "scf_algorithm": scf_algorithm,
        "pcm_dielectric": pcm_dielectric,
        "smd_solvent": smd_solvent,
        "overwrite_inputs": overwrite_inputs,
    }

    if scf_algorithm.lower() == "gdm" and "max_scf_cycles" not in qchem_input_params:
        qchem_input_params["max_scf_cycles"] = 200

    calc = QChem(
        input_atoms=atoms,
        cores=cores,
        charge=charge,
        spin_multiplicity=spin_multiplicity,
        qchem_input_params=qchem_input_params,
    )
    atoms.calc = calc
    atoms = run_calc(atoms)
    return summarize_run(
        atoms,
        input_atoms=input_atoms,
        charge_and_multiplicity=(checked_charge, checked_spin_multiplicity),
        additional_fields={"name": "Q-Chem Static"},
    )


@ct.electron
def relax_job(
    atoms: Atoms,
    cores: int | None = None,
    charge: int | None = None,
    spin_multiplicity: int | None = None,
    method: str = "wb97mv",
    basis: str = "def2-svpd",
    scf_algorithm: str = "diis",
    pcm_dielectric: str | None = None,
    smd_solvent: str | None = None,
    overwrite_inputs: dict | None = None,
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
    spin_multiplicity
        The spin multiplicity of the molecular system.
        Effectively defaults to the lowest spin state given the molecular structure and charge.
    method
        DFT exchange-correlation functional or other electronic structure method.
        Defaults to wB97M-V.
    basis
        Basis set.
        Defaults to def2-SVPD.
    scf_algorithm
        Algorithm used to converge the SCF.
        Defaults to "diis", but for particularly difficult cases, "gdm" should be employed instead.
    pcm_dielectric
        Dielectric constant of the optional polarizable continuum impicit solvation model.
        Defaults to None, in which case PCM will not be employed.
    smd_solvent
        Solvent to use for SMD implicit solvation model. Examples include "water", "ethanol", "methanol",
        and "acetonitrile". Refer to the Q-Chem manual for a complete list of solvents available.
        Defaults to None, in which case SMD will not be employed.
    overwrite_inputs
        Dictionary passed to pymatgen.io.qchem.QChemDictSet which can modify default values set therein
        as well as set additional Q-Chem parameters. See QChemDictSet documentation for more details.
    opt_swaps
        Dictionary of custom kwargs for run_ase_opt
            opt_defaults = {"fmax": 0.01, "max_steps": 1000, "optimizer": "Sella"}

    Returns
    -------
    dict
        Dictionary of results from quacc.schemas.ase.summarize_opt_run
    """

    # Reminder to self: exposing TRICs?

    checked_charge, checked_spin_multiplicity = check_charge_and_spin(atoms, charge, spin_multiplicity)

    opt_swaps = opt_swaps or {}
    opt_defaults = {
        "fmax": 0.01,
        "max_steps": 1000,
        "optimizer": Sella if Sella else FIRE,
        "optimizer_kwargs": {},
    }
    opt_flags = opt_defaults | opt_swaps
    if Sella:
        if opt_flags["optimizer"] == Sella and "order" not in opt_flags["optimizer_kwargs"]:
            opt_flags["optimizer_kwargs"]["order"] = 0

    if pcm_dielectric is not None and smd_solvent is not None:
        raise ValueError("PCM and SMD cannot be employed simultaneously! Exiting...")

    overwrite_inputs = overwrite_inputs or {}
    if "rem" not in overwrite_inputs:
        overwrite_inputs["rem"] = {}
    if "method" not in overwrite_inputs["rem"]:
        overwrite_inputs["rem"]["method"] = method

    qchem_input_params = {
        "basis_set": basis,
        "scf_algorithm": scf_algorithm,
        "pcm_dielectric": pcm_dielectric,
        "smd_solvent": smd_solvent,
        "overwrite_inputs": overwrite_inputs,
    }

    if scf_algorithm.lower() == "gdm" and "max_scf_cycles" not in qchem_input_params:
        qchem_input_params["max_scf_cycles"] = 200

    calc = QChem(
        input_atoms=atoms,
        cores=cores,
        charge=charge,
        spin_multiplicity=spin_multiplicity,
        qchem_input_params=qchem_input_params,
    )
    atoms.calc = calc
    dyn = run_ase_opt(atoms, **opt_flags)

    return summarize_opt_run(
        dyn,
        check_convergence=check_convergence,
        charge_and_multiplicity=(checked_charge, checked_spin_multiplicity),
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
    spin_multiplicity: int | None = None,
    method: str = "wb97mv",
    basis: str = "def2-svpd",
    scf_algorithm: str = "diis",
    pcm_dielectric: str | None = None,
    smd_solvent: str | None = None,
    overwrite_inputs: dict | None = None,
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
    spin_multiplicity
        The spin multiplicity of the molecular system.
        Effectively defaults to the lowest spin state given the molecular structure and charge.
    method
        DFT exchange-correlation functional or other electronic structure method.
        Defaults to wB97M-V.
    basis
        Basis set.
        Defaults to def2-SVPD.
    scf_algorithm
        Algorithm used to converge the SCF.
        Defaults to "diis", but for particularly difficult cases, "gdm" should be employed instead.
    pcm_dielectric
        Dielectric constant of the optional polarizable continuum impicit solvation model.
        Defaults to None, in which case PCM will not be employed.
    smd_solvent
        Solvent to use for SMD implicit solvation model. Examples include "water", "ethanol", "methanol",
        and "acetonitrile". Refer to the Q-Chem manual for a complete list of solvents available.
        Defaults to None, in which case SMD will not be employed.
    overwrite_inputs
        Dictionary passed to pymatgen.io.qchem.QChemDictSet which can modify default values set therein
        as well as set additional Q-Chem parameters. See QChemDictSet documentation for more details.
    opt_swaps
        Dictionary of custom kwargs for run_ase_opt
            opt_defaults = {"fmax": 0.01, "max_steps": 1000, "optimizer": "Sella"}

    Returns
    -------
    dict
        Dictionary of results from quacc.schemas.ase.summarize_opt_run
    """

    # Reminders to self:
    #   - exposing TRICs?
    #   - passing initial Hessian?

    checked_charge, checked_spin_multiplicity = check_charge_and_spin(atoms, charge, spin_multiplicity)

    opt_swaps = opt_swaps or {}
    opt_defaults = {
        "fmax": 0.01,
        "max_steps": 1000,
        "optimizer": Sella,
        "optimizer_kwargs": {},
    }
    opt_flags = opt_defaults | opt_swaps
    if opt_flags["optimizer"] != Sella:
        raise RuntimeError("Only Sella should be used for TS optimization! Exiting...")

    if pcm_dielectric is not None and smd_solvent is not None:
        raise ValueError("PCM and SMD cannot be employed simultaneously! Exiting...")

    overwrite_inputs = overwrite_inputs or {}
    if "rem" not in overwrite_inputs:
        overwrite_inputs["rem"] = {}
    if "method" not in overwrite_inputs["rem"]:
        overwrite_inputs["rem"]["method"] = method

    qchem_input_params = {
        "basis_set": basis,
        "scf_algorithm": scf_algorithm,
        "pcm_dielectric": pcm_dielectric,
        "smd_solvent": smd_solvent,
        "overwrite_inputs": overwrite_inputs,
    }

    if scf_algorithm.lower() == "gdm" and "max_scf_cycles" not in qchem_input_params:
        qchem_input_params["max_scf_cycles"] = 200

    calc = QChem(
        input_atoms=atoms,
        cores=cores,
        charge=charge,
        spin_multiplicity=spin_multiplicity,
        qchem_input_params=qchem_input_params,
    )
    atoms.calc = calc
    dyn = run_ase_opt(atoms, **opt_flags)

    return summarize_opt_run(
        dyn,
        check_convergence=check_convergence,
        charge_and_multiplicity=(checked_charge, checked_spin_multiplicity),
        additional_fields={"name": "Q-Chem TS Optimization"},
    )


@ct.electron
@requires(
    Sella,
    "Sella must be installed. pip install sella",
)
def irc_job(
    atoms: Atoms,
    direction: Literal["forward", "reverse"],
    cores: int | None = None,
    charge: int | None = None,
    spin_multiplicity: int | None = None,
    method: str = "wb97mv",
    basis: str = "def2-svpd",
    scf_algorithm: str = "diis",
    pcm_dielectric: str | None = None,
    smd_solvent: str | None = None,
    overwrite_inputs: dict | None = None,
    opt_swaps: dict | None = None,
    check_convergence: bool = True,
) -> dict:
    """
    IRC optimize a molecular structure.

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
    spin_multiplicity
        The spin multiplicity of the molecular system.
        Effectively defaults to the lowest spin state given the molecular structure and charge.
    method
        DFT exchange-correlation functional or other electronic structure method.
        Defaults to wB97M-V.
    basis
        Basis set.
        Defaults to def2-SVPD.
    scf_algorithm
        Algorithm used to converge the SCF.
        Defaults to "diis", but for particularly difficult cases, "gdm" should be employed instead.
    pcm_dielectric
        Dielectric constant of the optional polarizable continuum impicit solvation model.
        Defaults to None, in which case PCM will not be employed.
    smd_solvent
        Solvent to use for SMD implicit solvation model. Examples include "water", "ethanol", "methanol",
        and "acetonitrile". Refer to the Q-Chem manual for a complete list of solvents available.
        Defaults to None, in which case SMD will not be employed.
    overwrite_inputs
        Dictionary passed to pymatgen.io.qchem.QChemDictSet which can modify default values set therein
        as well as set additional Q-Chem parameters. See QChemDictSet documentation for more details.
    opt_swaps
        Dictionary of custom kwargs for run_ase_opt
            opt_defaults = {"fmax": 0.01, "max_steps": 1000, "optimizer": "Sella"}

    Returns
    -------
    dict
        Dictionary of results from quacc.schemas.ase.summarize_opt_run
    """

    # Reminders to self:
    #   - exposing TRICs?
    #   - passing initial Hessian?

    checked_charge, checked_spin_multiplicity = check_charge_and_spin(atoms, charge, spin_multiplicity)

    if direction not in ["forward", "reverse"]:
        raise ValueError("direction must be 'forward' or 'reverse'! Exiting...")

    opt_swaps = opt_swaps or {}
    opt_defaults = {
        "fmax": 0.01,
        "max_steps": 1000,
        "optimizer": IRC,
        "optimizer_kwargs": {},
        "run_kwargs": {"direction": direction},
    }
    opt_flags = opt_defaults | opt_swaps
    if opt_flags["optimizer"] != IRC:
        raise RuntimeError("Only Sella's IRC should be used for IRC optimization! Exiting...")

    if pcm_dielectric is not None and smd_solvent is not None:
        raise ValueError("PCM and SMD cannot be employed simultaneously! Exiting...")

    overwrite_inputs = overwrite_inputs or {}
    if "rem" not in overwrite_inputs:
        overwrite_inputs["rem"] = {}
    if "method" not in overwrite_inputs["rem"]:
        overwrite_inputs["rem"]["method"] = method

    qchem_input_params = {
        "basis_set": basis,
        "scf_algorithm": scf_algorithm,
        "pcm_dielectric": pcm_dielectric,
        "smd_solvent": smd_solvent,
        "overwrite_inputs": overwrite_inputs,
    }

    if scf_algorithm.lower() == "gdm" and "max_scf_cycles" not in qchem_input_params:
        qchem_input_params["max_scf_cycles"] = 200

    calc = QChem(
        input_atoms=atoms,
        cores=cores,
        charge=charge,
        spin_multiplicity=spin_multiplicity,
        qchem_input_params=qchem_input_params,
    )
    atoms.calc = calc
    dyn = run_ase_opt(atoms, **opt_flags)

    return summarize_opt_run(
        dyn,
        check_convergence=check_convergence,
        charge_and_multiplicity=(checked_charge, checked_spin_multiplicity),
        additional_fields={"name": "Q-Chem IRC Optimization"},
    )


@ct.electron
@requires(
    Sella,
    "Sella must be installed. pip install sella",
)
def quasi_irc_job(
    atoms: Atoms,
    direction: Literal["forward", "reverse"],
    cores: int | None = None,
    charge: int | None = None,
    spin_multiplicity: int | None = None,
    method: str = "wb97mv",
    basis: str = "def2-svpd",
    scf_algorithm: str = "diis",
    pcm_dielectric: str | None = None,
    smd_solvent: str | None = None,
    overwrite_inputs: dict | None = None,
    irc_swaps: dict | None = None,
    relax_swaps: dict | None = None,
    check_convergence: bool = True,
) -> dict:
    """
    Quasi-IRC optimize a molecular structure.

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
    spin_multiplicity
        The spin multiplicity of the molecular system.
        Effectively defaults to the lowest spin state given the molecular structure and charge.
    method
        DFT exchange-correlation functional or other electronic structure method.
        Defaults to wB97M-V.
    basis
        Basis set.
        Defaults to def2-SVPD.
    scf_algorithm
        Algorithm used to converge the SCF.
        Defaults to "diis", but for particularly difficult cases, "gdm" should be employed instead.
    pcm_dielectric
        Dielectric constant of the optional polarizable continuum impicit solvation model.
        Defaults to None, in which case PCM will not be employed.
    smd_solvent
        Solvent to use for SMD implicit solvation model. Examples include "water", "ethanol", "methanol",
        and "acetonitrile". Refer to the Q-Chem manual for a complete list of solvents available.
        Defaults to None, in which case SMD will not be employed.
    overwrite_inputs
        Dictionary passed to pymatgen.io.qchem.QChemDictSet which can modify default values set therein
        as well as set additional Q-Chem parameters. See QChemDictSet documentation for more details.
    irc_swaps
        Dictionary of custom kwargs for the irc_job.
    relax_swaps
        Dictionary of custom kwargs for the relax_job.

    Returns
    -------
    dict
        Dictionary of results from quacc.schemas.ase.summarize_opt_run
    """

    # Reminders to self:
    #   - exposing TRICs?
    #   - passing initial Hessian?

    if direction not in ["forward", "reverse"]:
        raise ValueError("direction must be 'forward' or 'reverse'! Exiting...")

    if pcm_dielectric is not None and smd_solvent is not None:
        raise ValueError("PCM and SMD cannot be employed simultaneously! Exiting...")

    irc_swaps = irc_swaps or {}
    irc_swaps_defaults = {
        "fmax": 100,
        "max_steps": 10,
    }
    irc_swaps = irc_swaps_defaults | irc_swaps

    irc_summary = irc_job(
        atoms=atoms,
        direction=direction,
        cores=cores,
        charge=charge,
        spin_multiplicity=spin_multiplicity,
        method=method,
        basis=basis,
        scf_algorithm=scf_algorithm,
        pcm_dielectric=pcm_dielectric,
        smd_solvent=smd_solvent,
        overwrite_inputs=overwrite_inputs,
        opt_swaps=irc_swaps,
        check_convergence=False,
    )

    relax_summary = relax_job(
        atoms=irc_summary["atoms"],
        cores=cores,
        charge=charge,
        spin_multiplicity=spin_multiplicity,
        method=method,
        basis=basis,
        scf_algorithm=scf_algorithm,
        pcm_dielectric=pcm_dielectric,
        smd_solvent=smd_solvent,
        overwrite_inputs=overwrite_inputs,
        opt_swaps=relax_swaps,
        check_convergence=check_convergence,
    )

    relax_summary["initial_irc"] = irc_summary

    return relax_summary

"""Core recipes for the Q-Chem"""
from __future__ import annotations

import multiprocessing
from typing import Literal

import covalent as ct
from ase.atoms import Atoms
from ase.optimize import FIRE
from monty.dev import requires

from quacc import SETTINGS
from quacc.calculators.qchem import QChem
from quacc.schemas.ase import OptSchema, RunSchema, summarize_opt_run, summarize_run
from quacc.util.atoms import check_charge_and_spin
from quacc.util.calc import run_ase_opt, run_calc
from quacc.util.dicts import remove_dict_empties

try:
    from sella import IRC, Sella

    has_sella = True
except ImportError:
    has_sella = False


@ct.electron
def static_job(
    atoms: Atoms | dict,
    charge: int | None = None,
    spin_multiplicity: int | None = None,
    method: str = "wb97mv",
    basis: str = "def2-tzvpd",
    scf_algorithm: str = "diis",
    pcm_dielectric: str | None = None,
    smd_solvent: str | None = None,
    overwrite_inputs: dict | None = None,
    n_cores: int | None = None,
) -> RunSchema:
    """
    Carry out a single-point calculation.

    Parameters
    ----------
    atoms
        Atoms object or a dictionary with the key "atoms" and an Atoms object as the value
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
    n_cores
        Number of cores to use for the Q-Chem calculation.
        Defaults to use all cores available on a given node.

    Returns
    -------
    RunSchema
        Dictionary of results from quacc.schemas.ase.summarize_run
    """
    atoms = atoms if isinstance(atoms, Atoms) else atoms["atoms"]
    checked_charge, checked_spin_multiplicity = check_charge_and_spin(
        atoms, charge, spin_multiplicity
    )

    qchem_defaults = {
        "method": method,
        "charge": charge,
        "spin_multiplicity": spin_multiplicity,
        "cores": n_cores or multiprocessing.cpu_count(),
        "qchem_input_params": {
            "basis_set": basis,
            "scf_algorithm": scf_algorithm,
            "pcm_dielectric": pcm_dielectric,
            "smd_solvent": smd_solvent,
            "overwrite_inputs": overwrite_inputs,
            "max_scf_cycles": 200 if scf_algorithm.lower() == "gdm" else None,
        },
    }
    qchem_flags = remove_dict_empties(qchem_defaults)

    atoms.calc = QChem(atoms, **qchem_flags)
    final_atoms = run_calc(atoms)

    return summarize_run(
        final_atoms,
        input_atoms=atoms,
        charge_and_multiplicity=(checked_charge, checked_spin_multiplicity),
        additional_fields={"name": "Q-Chem Static"},
    )


@ct.electron
def relax_job(
    atoms: Atoms | dict,
    charge: int | None = None,
    spin_multiplicity: int | None = None,
    method: str = "wb97mv",
    basis: str = "def2-svpd",
    scf_algorithm: str = "diis",
    pcm_dielectric: str | None = None,
    smd_solvent: str | None = None,
    overwrite_inputs: dict | None = None,
    opt_swaps: dict | None = None,
    n_cores: int | None = None,
) -> OptSchema:
    """
    Optimize aka "relax" a molecular structure.

    Parameters
    ----------
    atoms
        Atoms object or a dictionary with the key "atoms" and an Atoms object as the value
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
    n_cores
        Number of cores to use for the Q-Chem calculation.
        Defaults to use all cores available on a given node.

    Returns
    -------
    OptSchema
        Dictionary of results from quacc.schemas.ase.summarize_opt_run
    """

    # TODO: exposing TRICs?
    atoms = atoms if isinstance(atoms, Atoms) else atoms["atoms"]
    checked_charge, checked_spin_multiplicity = check_charge_and_spin(
        atoms, charge, spin_multiplicity
    )

    qchem_defaults = {
        "method": method,
        "charge": charge,
        "spin_multiplicity": spin_multiplicity,
        "cores": n_cores or multiprocessing.cpu_count(),
        "qchem_input_params": {
            "basis_set": basis,
            "scf_algorithm": scf_algorithm,
            "pcm_dielectric": pcm_dielectric,
            "smd_solvent": smd_solvent,
            "overwrite_inputs": overwrite_inputs,
            "max_scf_cycles": 200 if scf_algorithm.lower() == "gdm" else None,
        },
    }
    qchem_flags = remove_dict_empties(qchem_defaults)

    opt_swaps = opt_swaps or {}
    opt_defaults = {
        "fmax": 0.01,
        "max_steps": 1000,
        "optimizer": Sella if has_sella else FIRE,
        "optimizer_kwargs": None,
    }
    if opt_flags["optimizer"] == Sella and "order" not in opt_flags["optimizer_kwargs"]:
        opt_flags["optimizer_kwargs"]["order"] = 0
    opt_flags = remove_dict_empties(opt_defaults | opt_swaps)

    atoms.calc = QChem(atoms, **qchem_flags)
    dyn = run_ase_opt(atoms, **opt_flags)

    return summarize_opt_run(
        dyn,
        charge_and_multiplicity=(checked_charge, checked_spin_multiplicity),
        additional_fields={"name": "Q-Chem Optimization"},
    )


@ct.electron
@requires(
    has_sella,
    "Sella must be installed. Try pip install sella.",
)
def ts_job(
    atoms: Atoms | dict,
    charge: int | None = None,
    spin_multiplicity: int | None = None,
    method: str = "wb97mv",
    basis: str = "def2-svpd",
    scf_algorithm: str = "diis",
    pcm_dielectric: str | None = None,
    smd_solvent: str | None = None,
    overwrite_inputs: dict | None = None,
    opt_swaps: dict | None = None,
    n_cores: int | None = None,
) -> OptSchema:
    """
    TS optimize a molecular structure.

    Parameters
    ----------
    atoms
        Atoms object or a dictionary with the key "atoms" and an Atoms object as the value
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
    n_cores
        Number of cores to use for the Q-Chem calculation.
        Defaults to use all cores available on a given node.

    Returns
    -------
    OptSchema
        Dictionary of results from quacc.schemas.ase.summarize_opt_run
    """

    # TODO:
    #   - exposing TRICs?
    #   - passing initial Hessian?
    atoms = atoms if isinstance(atoms, Atoms) else atoms["atoms"]
    checked_charge, checked_spin_multiplicity = check_charge_and_spin(
        atoms, charge, spin_multiplicity
    )

    qchem_defaults = {
        "method": method,
        "charge": charge,
        "spin_multiplicity": spin_multiplicity,
        "cores": n_cores or multiprocessing.cpu_count(),
        "qchem_input_params": {
            "basis_set": basis,
            "scf_algorithm": scf_algorithm,
            "pcm_dielectric": pcm_dielectric,
            "smd_solvent": smd_solvent,
            "overwrite_inputs": overwrite_inputs,
            "max_scf_cycles": 200 if scf_algorithm.lower() == "gdm" else None,
        },
    }
    qchem_flags = remove_dict_empties(qchem_defaults)

    opt_swaps = opt_swaps or {}
    opt_defaults = {
        "fmax": 0.01,
        "max_steps": 1000,
        "optimizer": Sella,
        "optimizer_kwargs": None,
    }
    opt_flags = remove_dict_empties(opt_defaults | opt_swaps)
    if opt_flags["optimizer"] != Sella:
        raise ValueError("Only Sella should be used for TS optimization.")

    atoms.calc = QChem(atoms, **qchem_flags)
    dyn = run_ase_opt(atoms, **opt_flags)

    return summarize_opt_run(
        dyn,
        charge_and_multiplicity=(checked_charge, checked_spin_multiplicity),
        additional_fields={"name": "Q-Chem TS Optimization"},
    )


@ct.electron
@requires(
    has_sella,
    "Sella must be installed. Try pip install sella.",
)
def irc_job(
    atoms: Atoms | dict,
    direction: Literal["forward", "reverse"],
    charge: int | None = None,
    spin_multiplicity: int | None = None,
    method: str = "wb97mv",
    basis: str = "def2-svpd",
    scf_algorithm: str = "diis",
    pcm_dielectric: str | None = None,
    smd_solvent: str | None = None,
    overwrite_inputs: dict | None = None,
    opt_swaps: dict | None = None,
    n_cores: int | None = None,
) -> OptSchema:
    """
    IRC optimize a molecular structure.

    Parameters
    ----------
    atoms
        Atoms object or a dictionary with the key "atoms" and an Atoms object as the value
    direction
        Direction of the IRC. Should be "forward" or "reverse".
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
    n_cores
        Number of cores to use for the Q-Chem calculation.
        Defaults to use all cores available on a given node.

    Returns
    -------
    OptSchema
        Dictionary of results from quacc.schemas.ase.summarize_opt_run
    """

    # TODO: 1) expose TRICs?; 2) passing initial Hessian?
    atoms = atoms if isinstance(atoms, Atoms) else atoms["atoms"]
    checked_charge, checked_spin_multiplicity = check_charge_and_spin(
        atoms, charge, spin_multiplicity
    )

    qchem_defaults = {
        "method": method,
        "charge": charge,
        "spin_multiplicity": spin_multiplicity,
        "cores": n_cores or multiprocessing.cpu_count(),
        "qchem_input_params": {
            "basis_set": basis,
            "scf_algorithm": scf_algorithm,
            "pcm_dielectric": pcm_dielectric,
            "smd_solvent": smd_solvent,
            "overwrite_inputs": overwrite_inputs,
            "max_scf_cycles": 200 if scf_algorithm.lower() == "gdm" else None,
        },
    }
    qchem_flags = remove_dict_empties(qchem_defaults)

    opt_swaps = opt_swaps or {}
    opt_defaults = {
        "fmax": 0.01,
        "max_steps": 1000,
        "optimizer": IRC,
        "optimizer_kwargs": {"keep_going": True},
        "run_kwargs": {"direction": direction},
    }
    opt_flags = remove_dict_empties(opt_defaults | opt_swaps)
    if opt_flags["optimizer"] != IRC:
        raise ValueError("Only Sella's IRC should be used for IRC optimization.")

    atoms.calc = QChem(atoms, **qchem_flags)
    dyn = run_ase_opt(atoms, **opt_flags)

    return summarize_opt_run(
        dyn,
        charge_and_multiplicity=(checked_charge, checked_spin_multiplicity),
        additional_fields={"name": "Q-Chem IRC Optimization"},
    )


@ct.electron
@requires(
    Sella,
    "Sella must be installed. Try pip install sella.",
)
def quasi_irc_job(
    atoms: Atoms | dict,
    direction: Literal["forward", "reverse"],
    common_kwargs: dict | None = None,
    irc_opt_swaps: dict | None = None,
    relax_opt_swaps: dict | None = None,
) -> OptSchema:
    """
    Quasi-IRC optimize a molecular structure.

    Parameters
    ----------
    atoms
        Atoms object.
    direction
        Direction of the IRC. Should be "forward" or "reverse".
    common_kwargs
        Dictionary of kwargs that are passed as input to both irc_job and relax_job.
    irc_opt_swaps
        Dictionary of opt_swap kwargs for the irc_job.
    relax_opt_swaps
        Dictionary of opt_swap kwargs for the relax_job.

    Returns
    -------
    OptSchema
        Dictionary of results from quacc.schemas.ase.summarize_opt_run
    """

    common_kwargs = common_kwargs or {}
    irc_opt_swaps = irc_opt_swaps or {}
    relax_opt_swaps = relax_opt_swaps or {}
    default_convergence = SETTINGS.CHECK_CONVERGENCE

    irc_opt_swaps_defaults = {
        "fmax": 100,
        "max_steps": 10,
    }
    irc_opt_swaps = remove_dict_empties(irc_opt_swaps_defaults | irc_opt_swaps)

    SETTINGS.CHECK_CONVERGENCE = False
    irc_summary = irc_job(
        atoms,
        direction=direction,
        opt_swaps=irc_opt_swaps,
        **common_kwargs,
    )

    SETTINGS.CHECK_CONVERGENCE = default_convergence
    relax_summary = relax_job(
        irc_summary,
        opt_swaps=relax_opt_swaps,
        **common_kwargs,
    )

    relax_summary["initial_irc"] = irc_summary

    return relax_summary

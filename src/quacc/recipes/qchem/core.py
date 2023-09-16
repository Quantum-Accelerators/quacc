"""Core recipes for the Q-Chem"""
from __future__ import annotations

import multiprocessing
from typing import TYPE_CHECKING

from ase.optimize import FIRE

from quacc import job
from quacc.calculators.qchem import QChem
from quacc.schemas import fetch_atoms
from quacc.schemas.ase import summarize_opt_run
from quacc.schemas.qchem import QchemSchema, summarize_run
from quacc.utils.calc import run_ase_opt, run_calc
from quacc.utils.dicts import merge_dicts, remove_dict_empties

try:
    from sella import Sella

    has_sella = True
except ImportError:
    has_sella = False

if TYPE_CHECKING:
    from ase import Atoms

    from quacc.schemas.ase import OptSchema, RunSchema


@job
def static_job(
    atoms: Atoms | dict,
    charge: int,
    spin_multiplicity: int,
    method: str = "wb97mv",
    basis: str = "def2-tzvpd",
    scf_algorithm: str = "diis",
    pcm_dielectric: str | None = None,
    smd_solvent: str | None = None,
    n_cores: int | None = None,
    overwrite_inputs: dict | None = None,
    copy_files: list[str] | None = None,
) -> QchemSchema:
    """
    Carry out a single-point calculation.

    Parameters
    ----------
    atoms
        Atoms object or a dictionary with the key "atoms" and an Atoms object as
        the value
    charge
        Charge of the system.
    spin_multiplicity
        Multiplicity of the system.
    method
        DFT exchange-correlation functional or other electronic structure
        method. Defaults to wB97M-V.
    basis
        Basis set. Defaults to def2-TZVPD.
    scf_algorithm
        Algorithm used to converge the SCF. Defaults to "diis", but for
        particularly difficult cases, "gdm" should be employed instead.
    pcm_dielectric
        Dielectric constant of the optional polarizable continuum impicit
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
    copy_files
        Files to copy to the runtime directory.

    Returns
    -------
    QchemSchema
        Dictionary of results from [quacc.schemas.qchem.summarize_run][]
    """
    atoms = fetch_atoms(atoms)

    qchem_defaults = {
        "basis_set": basis,
        "scf_algorithm": scf_algorithm,
        "method": method,
        "charge": charge,
        "spin_multiplicity": spin_multiplicity,
        "cores": n_cores or multiprocessing.cpu_count(),
        "qchem_input_params": {
            "pcm_dielectric": pcm_dielectric,
            "smd_solvent": smd_solvent,
            "overwrite_inputs": overwrite_inputs,
            "max_scf_cycles": 200 if scf_algorithm.lower() == "gdm" else None,
        },
    }
    qchem_flags = remove_dict_empties(qchem_defaults)

    atoms.calc = QChem(atoms, **qchem_flags)
    atoms = run_calc(atoms, copy_files=copy_files)

    return summarize_run(atoms, additional_fields={"name": "Q-Chem Static"})


@job
def relax_job(
    atoms: Atoms | dict,
    charge: int,
    spin_multiplicity: int,
    method: str = "wb97mv",
    basis: str = "def2-svpd",
    scf_algorithm: str = "diis",
    pcm_dielectric: str | None = None,
    smd_solvent: str | None = None,
    n_cores: int | None = None,
    overwrite_inputs: dict | None = None,
    opt_swaps: dict | None = None,
    copy_files: list[str] | None = None,
) -> OptSchema:
    """
    Optimize aka "relax" a molecular structure with an ASE optimizer.

    Parameters
    ----------
    atoms
        Atoms object or a dictionary with the key "atoms" and an Atoms object as
        the value
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
        Dielectric constant of the optional polarizable continuum impicit
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
    opt_swaps
        Dictionary of custom kwargs for [quacc.utils.calc.run_ase_opt][]

        ???+ Note

             Overrides the following defaults:

            ```python
            {"fmax": 0.01, "max_steps": 1000, "optimizer": Sella}
            ```
    copy_files
        Files to copy to the runtime directory.

    Returns
    -------
    OptSchema
        Dictionary of results from [quacc.schemas.ase.summarize_opt_run][]
    """

    # TODO: exposing TRICs?
    atoms = fetch_atoms(atoms)

    qchem_defaults = {
        "basis_set": basis,
        "scf_algorithm": scf_algorithm,
        "method": method,
        "charge": charge,
        "spin_multiplicity": spin_multiplicity,
        "cores": n_cores or multiprocessing.cpu_count(),
        "qchem_input_params": {
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
    }
    opt_flags = merge_dicts(opt_defaults, opt_swaps)
    if opt_flags["optimizer"].__name__ == "Sella" and "order" not in opt_flags.get(
        "optimizer_kwargs", {}
    ):
        if "optimizer_kwargs" not in opt_flags:
            opt_flags["optimizer_kwargs"] = {}
        opt_flags["optimizer_kwargs"]["order"] = 0

    atoms.calc = QChem(atoms, **qchem_flags)
    dyn = run_ase_opt(atoms, copy_files=copy_files, **opt_flags)

    return summarize_opt_run(
        dyn,
        charge_and_multiplicity=(charge, spin_multiplicity),
        additional_fields={"name": "Q-Chem Optimization"},
    )


@job
def internal_relax_job(
    atoms: Atoms | dict,
    charge: int,
    spin_multiplicity: int,
    method: str = "wb97mv",
    basis: str = "def2-svpd",
    scf_algorithm: str = "diis",
    pcm_dielectric: str | None = None,
    smd_solvent: str | None = None,
    n_cores: int | None = None,
    overwrite_inputs: dict | None = None,
    copy_files: list[str] | None = None,
) -> RunSchema:
    """
    Optimize aka "relax" a molecular structure with Q-Chem optimizers.

    Parameters
    ----------
    atoms
        Atoms object or a dictionary with the key "atoms" and an Atoms object as
        the value
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
        Dielectric constant of the optional polarizable continuum impicit
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
    copy_files
        Files to copy to the runtime directory.

    Returns
    -------
    RunSchema
        Dictionary of results from [quacc.schemas.ase.summarize_run][]
    """
    atoms = fetch_atoms(atoms)

    qchem_defaults = {
        "job_type": "opt",
        "basis_set": basis,
        "scf_algorithm": scf_algorithm,
        "method": method,
        "charge": charge,
        "spin_multiplicity": spin_multiplicity,
        "cores": n_cores or multiprocessing.cpu_count(),
        "qchem_input_params": {
            "pcm_dielectric": pcm_dielectric,
            "smd_solvent": smd_solvent,
            "overwrite_inputs": overwrite_inputs,
            "max_scf_cycles": 200 if scf_algorithm.lower() == "gdm" else None,
        },
    }
    qchem_flags = remove_dict_empties(qchem_defaults)

    atoms.calc = QChem(atoms, **qchem_flags)
    final_atoms = run_calc(atoms, copy_files=copy_files)

    return summarize_run(
        final_atoms,
        input_atoms=atoms,
        charge_and_multiplicity=(charge, spin_multiplicity),
        additional_fields={"name": "Q-Chem Optimization (Internal)"},
    )


@job
def freq_job(
    atoms: Atoms | dict,
    charge: int,
    spin_multiplicity: int,
    method: str = "wb97mv",
    basis: str = "def2-svpd",
    scf_algorithm: str = "diis",
    pcm_dielectric: str | None = None,
    smd_solvent: str | None = None,
    n_cores: int | None = None,
    overwrite_inputs: dict | None = None,
    copy_files: list[str] | None = None,
) -> QchemSchema:
    """
    Perform a frequency calculation on a molecular structure.

    Parameters
    ----------
    atoms
        Atoms object or a dictionary with the key "atoms" and an Atoms object as
        the value
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
        Dielectric constant of the optional polarizable continuum impicit
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
    copy_files
        Files to copy to the runtime directory.

    Returns
    -------
    RunSchema
        Dictionary of results from [quacc.schemas.ase.summarize_run][]
    """

    atoms = fetch_atoms(atoms)

    qchem_defaults = {
        "job_type": "freq",
        "basis_set": basis,
        "scf_algorithm": scf_algorithm,
        "method": method,
        "charge": charge,
        "spin_multiplicity": spin_multiplicity,
        "cores": n_cores or multiprocessing.cpu_count(),
        "qchem_input_params": {
            "pcm_dielectric": pcm_dielectric,
            "smd_solvent": smd_solvent,
            "overwrite_inputs": overwrite_inputs,
            "max_scf_cycles": 200 if scf_algorithm.lower() == "gdm" else None,
        },
    }
    qchem_flags = remove_dict_empties(qchem_defaults)

    atoms.calc = QChem(atoms, **qchem_flags)
    atoms = run_calc(atoms, copy_files=copy_files)

    return summarize_run(atoms, additional_fields={"name": "Q-Chem Frequency"})

"""Core recipes for the Q-Chem."""
from __future__ import annotations

import multiprocessing
from typing import TYPE_CHECKING

from ase.optimize import FIRE

from quacc import SETTINGS, job
from quacc.calculators.qchem import QChem
from quacc.runners.calc import run_ase_calc, run_ase_opt
from quacc.schemas.ase import summarize_opt_run, summarize_run
from quacc.utils.dicts import merge_dicts, remove_dict_nones

try:
    from sella import Sella

    has_sella = True
except ImportError:
    has_sella = False

if TYPE_CHECKING:
    from typing import Any

    from ase import Atoms

    from quacc.schemas.ase import OptSchema, RunSchema


@job
def static_job(
    atoms: Atoms,
    charge: int,
    spin_multiplicity: int,
    method: str = "wb97mv",
    basis: str = "def2-tzvpd",
    scf_algorithm: str = "diis",
    pcm_dielectric: str | None = None,
    smd_solvent: str | None = None,
    n_cores: int | None = None,
    overwrite_inputs: dict[str, Any] | None = None,
    copy_files: list[str] | None = None,
) -> RunSchema:
    """
    Carry out a single-point calculation.

    ??? Note

        Calculator Defaults:

        ```python
        {
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
            "nbo_params": {"version": 7} if SETTINGS.QCHEM_NBO_EXE else None,
            },
        }
        ```

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
    RunSchema
        Dictionary of results from [quacc.schemas.ase.summarize_run][]
    """
    defaults = {
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
            "nbo_params": {"version": 7} if SETTINGS.QCHEM_NBO_EXE else None,
        },
    }

    return _base_job(
        atoms,
        charge,
        spin_multiplicity,
        defaults=defaults,
        additional_fields={"name": "Q-Chem Static"},
        copy_files=copy_files,
    )


@job
def internal_relax_job(
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
    copy_files: list[str] | None = None,
) -> RunSchema:
    """
    Optimize aka "relax" a molecular structure with Q-Chem optimizers.

    ??? Note

        Calculator Defaults:

        ```python
        {
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
        ```

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

    defaults = {
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
    return _base_job(
        atoms,
        charge,
        spin_multiplicity,
        defaults=defaults,
        additional_fields={"name": "Q-Chem Optimization (Internal)"},
        copy_files=copy_files,
    )


@job
def freq_job(
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
    copy_files: list[str] | None = None,
) -> RunSchema:
    """
    Perform a frequency calculation on a molecular structure.

    ??? Note

        Calculator Defaults:

        ```python
        {
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
        ```

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

    defaults = {
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
    return _base_job(
        atoms,
        charge,
        spin_multiplicity,
        defaults=defaults,
        copy_files=copy_files,
        additional_fields={"name": "Q-Chem Frequency"},
    )


@job
def relax_job(
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
    opt_swaps: dict[str, Any] | None = None,
    copy_files: list[str] | None = None,
) -> OptSchema:
    """
    Optimize aka "relax" a molecular structure with an ASE optimizer.

    ??? Note

        Calculator defaults:

        ```python
        {
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
        ```

        Optimizer defaults:

        ```python
        {"fmax": 0.01, "max_steps": 1000, "optimizer": Sella if has_sella else FIRE, "optimizer_kwargs": {"use_TRICs": False}}
        ```

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
        Dictionary of custom kwargs for [quacc.runners.calc.run_ase_opt][]
    copy_files
        Files to copy to the runtime directory.

    Returns
    -------
    OptSchema
        Dictionary of results from [quacc.schemas.ase.summarize_opt_run][]
    """

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
    opt_defaults = {
        "fmax": 0.01,
        "max_steps": 1000,
        "optimizer": Sella if has_sella else FIRE,
        "optimizer_kwargs": {"use_TRICs": False},
    }

    return _base_opt_job(
        atoms,
        charge,
        spin_multiplicity,
        qchem_defaults=qchem_defaults,
        opt_defaults=opt_defaults,
        opt_swaps=opt_swaps,
        additional_fields={"name": "Q-Chem Optimization"},
        copy_files=copy_files,
    )


def _base_job(
    atoms: Atoms,
    charge: int,
    spin_multiplicity: int,
    defaults: dict[str, Any] | None = None,
    additional_fields: dict[str, Any] | None = None,
    copy_files: list[str] | None = None,
) -> RunSchema:
    """
    Base job function used for Q-Chem recipes that don't rely on ASE optimizers or other
    ASE dynamics classes.

    Parameters
    ----------
    atoms
        Atoms object
    charge
        Charge of the system.
    spin_multiplicity
        Multiplicity of the system.
    defaults
        The default parameters for the recipe.
    additional_fields
        Any additional fields to set in the summary.
    copy_files
        Files to copy to the runtime directory.

    Returns
    -------
    RunSchema
        Dictionary of results from [quacc.schemas.ase.summarize_run][]
    """

    qchem_flags = remove_dict_nones(defaults)

    atoms.calc = QChem(atoms, **qchem_flags)
    final_atoms = run_ase_calc(atoms, copy_files=copy_files)

    return summarize_run(
        final_atoms,
        input_atoms=atoms,
        charge_and_multiplicity=(charge, spin_multiplicity),
        additional_fields=additional_fields,
    )


def _base_opt_job(
    atoms: Atoms,
    charge: int,
    spin_multiplicity: int,
    qchem_defaults: dict[str, Any] | None = None,
    opt_defaults: dict[str, Any] | None = None,
    opt_swaps: dict[str, Any] | None = None,
    additional_fields: dict[str, Any] | None = None,
    copy_files: list[str] | None = None,
) -> OptSchema:
    """
    Base function for Q-Chem recipes that involve ASE optimizers.

    Parameters
    ----------
    atoms
        Atoms object
    charge
        Charge of the system.
    spin_multiplicity
        Multiplicity of the system.
    qchem_defaults
        Default arguments for the Q-Chem calculator.
    opt_defaults
        Default arguments for the ASE optimizer.
    opt_swaps
        Dictionary of custom kwargs for [quacc.runners.calc.run_ase_opt][]
    copy_files
        Files to copy to the runtime directory.

    Returns
    -------
    OptSchema
        Dictionary of results from [quacc.schemas.ase.summarize_opt_run][]
    """
    # TODO:
    #   - passing initial Hessian?

    qchem_flags = remove_dict_nones(qchem_defaults)
    opt_flags = merge_dicts(opt_defaults, opt_swaps)

    atoms.calc = QChem(atoms, **qchem_flags)
    dyn = run_ase_opt(atoms, copy_files=copy_files, **opt_flags)

    return summarize_opt_run(
        dyn,
        charge_and_multiplicity=(charge, spin_multiplicity),
        additional_fields=additional_fields,
    )

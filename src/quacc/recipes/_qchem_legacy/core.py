"""Core recipes for the Q-Chem."""
from __future__ import annotations

from typing import TYPE_CHECKING

import psutil
from ase.optimize import FIRE

from quacc import SETTINGS, job
from quacc.recipes._qchem_legacy._base import base_fn, base_opt_fn

try:
    from sella import Sella

    has_sella = True
except ImportError:
    has_sella = False

if TYPE_CHECKING:
    from pathlib import Path
    from typing import Any

    from ase.atoms import Atoms

    from quacc.schemas._aliases.ase import OptSchema, RunSchema


@job
def static_job(
    atoms: Atoms,
    charge: int = 0,
    spin_multiplicity: int = 1,
    method: str = "wb97mv",
    basis: str = "def2-tzvpd",
    scf_algorithm: str = "diis",
    pcm_dielectric: str | None = None,
    smd_solvent: str | None = None,
    n_cores: int | None = None,
    overwrite_inputs: dict[str, Any] | None = None,
    copy_files: str | Path | list[str | Path] | None = None,
) -> RunSchema:
    """
    Carry out a single-point calculation.

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
    copy_files
        File(s) to copy to the runtime directory. If a directory is provided, it will be recursively unpacked.

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
        "cores": n_cores or psutil.cpu_count(logical=False),
        "qchem_input_params": {
            "pcm_dielectric": pcm_dielectric,
            "smd_solvent": smd_solvent,
            "overwrite_inputs": overwrite_inputs,
            "max_scf_cycles": 200 if scf_algorithm.lower() == "gdm" else 100,
            "nbo_params": {"version": 7} if SETTINGS.QCHEM_NBO_EXE else None,
        },
    }

    return base_fn(
        atoms,
        charge,
        spin_multiplicity,
        defaults=defaults,
        additional_fields={"name": "Q-Chem Static"},
        copy_files=copy_files,
    )


@job
def freq_job(
    atoms: Atoms,
    charge: int = 0,
    spin_multiplicity: int = 1,
    method: str = "wb97mv",
    basis: str = "def2-svpd",
    scf_algorithm: str = "diis",
    pcm_dielectric: str | None = None,
    smd_solvent: str | None = None,
    n_cores: int | None = None,
    overwrite_inputs: dict[str, Any] | None = None,
    copy_files: str | Path | list[str | Path] | None = None,
) -> RunSchema:
    """
    Perform a frequency calculation on a molecular structure.

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
    copy_files
        File(s) to copy to the runtime directory. If a directory is provided, it will be recursively unpacked.

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
        "cores": n_cores or psutil.cpu_count(logical=False),
        "qchem_input_params": {
            "pcm_dielectric": pcm_dielectric,
            "smd_solvent": smd_solvent,
            "overwrite_inputs": overwrite_inputs,
            "max_scf_cycles": 200 if scf_algorithm.lower() == "gdm" else 100,
        },
    }
    return base_fn(
        atoms,
        charge=charge,
        spin_multiplicity=spin_multiplicity,
        defaults=defaults,
        copy_files=copy_files,
        additional_fields={"name": "Q-Chem Frequency"},
    )


@job
def relax_job(
    atoms: Atoms,
    charge: int = 0,
    spin_multiplicity: int = 1,
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
    Optimize aka "relax" a molecular structure with an ASE optimizer.

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
        "optimizer": Sella if has_sella else FIRE,
    }

    return base_opt_fn(
        atoms,
        charge=charge,
        spin_multiplicity=spin_multiplicity,
        calc_defaults=calc_defaults,
        opt_defaults=opt_defaults,
        opt_params=opt_params,
        additional_fields={"name": "Q-Chem Optimization"},
        copy_files=copy_files,
    )

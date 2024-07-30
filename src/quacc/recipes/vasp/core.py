"""Core recipes for VASP."""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING, Literal
from warnings import warn

import numpy as np
from monty.os.path import zpath
from pymatgen.io.vasp import Vasprun

from quacc import flow, job
from quacc.recipes.vasp._base import (
    run_and_summarize,
    run_and_summarize_opt,
    run_and_summarize_vib_and_thermo,
)

if TYPE_CHECKING:
    from typing import Any

    from ase.atoms import Atoms

    from quacc.types import (
        DoubleRelaxSchema,
        Filenames,
        OptParams,
        SourceDirectory,
        VaspASEOptSchema,
        VaspSchema,
        VibKwargs,
        VibThermoSchema,
    )


@job
def static_job(
    atoms: Atoms,
    preset: str | None = "BulkSet",
    copy_files: SourceDirectory | dict[SourceDirectory, Filenames] | None = None,
    **calc_kwargs,
) -> VaspSchema:
    """
    Carry out a single-point calculation.

    Parameters
    ----------
    atoms
        Atoms object
    preset
        Preset to use from `quacc.calculators.vasp.presets`.
    copy_files
        Files to copy (and decompress) from source to the runtime directory.
    **calc_kwargs
        Custom kwargs for the Vasp calculator. Set a value to
        `None` to remove a pre-existing key entirely. For a list of available
        keys, refer to [quacc.calculators.vasp.vasp.Vasp][].

    Returns
    -------
    VaspSchema
        Dictionary of results from [quacc.schemas.vasp.VaspSummarize.run][].
        See the type-hint for the data structure.
    """
    calc_defaults = {
        "ismear": -5,
        "laechg": True,
        "lcharg": True,
        "lreal": False,
        "lwave": True,
        "nedos": 3001,
        "nsw": 0,
    }
    return run_and_summarize(
        atoms,
        preset=preset,
        calc_defaults=calc_defaults,
        calc_swaps=calc_kwargs,
        additional_fields={"name": "VASP Static"},
        copy_files=copy_files,
    )


@job
def relax_job(
    atoms: Atoms,
    preset: str | None = "BulkSet",
    relax_cell: bool = True,
    copy_files: SourceDirectory | dict[SourceDirectory, Filenames] | None = None,
    **calc_kwargs,
) -> VaspSchema:
    """
    Relax a structure.

    Parameters
    ----------
    atoms
        Atoms object
    preset
        Preset to use from `quacc.calculators.vasp.presets`.
    relax_cell
        True if a volume relaxation (ISIF = 3) should be performed. False if
        only the positions (ISIF = 2) should be updated.
    copy_files
        Files to copy (and decompress) from source to the runtime directory.
    **calc_kwargs
        Custom kwargs for the Vasp calculator. Set a value to
        `None` to remove a pre-existing key entirely. For a list of available
        keys, refer to the [quacc.calculators.vasp.vasp.Vasp][] calculator.

    Returns
    -------
    VaspSchema
        Dictionary of results from [quacc.schemas.vasp.VaspSummarize.run][].
        See the type-hint for the data structure.
    """
    if relax_cell:
        warn(
            "The `relax_cell` parameter will default to `False` by default in a future version for internal consistency throughout quacc. Please set `relax_cell=True` directly.",
            DeprecationWarning,
            stacklevel=3,
        )
    calc_defaults = {
        "ediffg": -0.02,
        "isif": 3 if relax_cell else 2,
        "ibrion": 2,
        "isym": 0,
        "lcharg": False,
        "lwave": False,
        "nsw": 200,
        "symprec": 1e-8,
    }
    return run_and_summarize(
        atoms,
        preset=preset,
        calc_defaults=calc_defaults,
        calc_swaps=calc_kwargs,
        additional_fields={"name": "VASP Relax"},
        copy_files=copy_files,
    )


@flow
def double_relax_flow(
    atoms: Atoms,
    preset: str | None = "BulkSet",
    relax_cell: bool = True,
    relax1_kwargs: dict[str, Any] | None = None,
    relax2_kwargs: dict[str, Any] | None = None,
) -> DoubleRelaxSchema:
    """
    Double-relax a structure. This is particularly useful for a few reasons:

    1. To carry out a cheaper pre-relaxation before the high-quality run.

    2. To carry out a GGA calculation before a meta-GGA or hybrid calculation
    that requires the GGA wavefunction.

    3. To carry out volume relaxations where large changes in volume
    can require a second relaxation to resolve forces.

    Parameters
    ----------
    atoms
        Atoms object
    preset
        Preset to use from `quacc.calculators.vasp.presets`.
    relax_cell
        True if a volume relaxation (ISIF = 3) should be performed. False if
        only the positions (ISIF = 2) should be updated.
    relax1_kwargs
        Dictionary of custom kwargs for the first relaxation.
    relax2_kwargs
        Dictionary of custom kwargs for the second relaxation.

    Returns
    -------
    DoubleRelaxSchema
        Dictionary of results from each step.
    """
    relax1_kwargs = relax1_kwargs or {}
    relax2_kwargs = relax2_kwargs or {}

    # Run first relaxation
    summary1 = relax_job(atoms, preset=preset, relax_cell=relax_cell, **relax1_kwargs)

    # Run second relaxation
    summary2 = relax_job(
        summary1["atoms"],
        preset=preset,
        relax_cell=relax_cell,
        copy_files={summary1["dir_name"]: ["WAVECAR*"]},
        **relax2_kwargs,
    )

    return {"relax1": summary1, "relax2": summary2}


@job
def ase_relax_job(
    atoms: Atoms,
    preset: str | None = "BulkSet",
    relax_cell: bool = True,
    opt_params: OptParams | None = None,
    copy_files: SourceDirectory | dict[SourceDirectory, Filenames] | None = None,
    **calc_kwargs,
) -> VaspASEOptSchema:
    """
    Relax a structure.

    Parameters
    ----------
    atoms
        Atoms object
    preset
        Preset to use from `quacc.calculators.vasp.presets`.
    relax_cell
        True if a volume relaxation should be performed. False if only the positions
        should be updated.
    opt_params
        Dictionary of custom kwargs for the optimization process. For a list
        of available keys, refer to [quacc.runners.ase.Runner.run_opt][].
    copy_files
        Files to copy (and decompress) from source to the runtime directory.
    **calc_kwargs
        Custom kwargs for the Vasp calculator. Set a value to
        `None` to remove a pre-existing key entirely. For a list of available
        keys, refer to the [quacc.calculators.vasp.vasp.Vasp][] calculator.

    Returns
    -------
    VaspASEOptSchema
        Dictionary of results. See the type-hint for the data structure.
    """
    if relax_cell:
        warn(
            "The `relax_cell` parameter will default to `False` by default in a future version for internal consistency throughout quacc. Please set `relax_cell=True` directly.",
            DeprecationWarning,
            stacklevel=3,
        )
    calc_defaults = {"lcharg": False, "lwave": False, "nsw": 0}
    opt_defaults = {"relax_cell": relax_cell}
    return run_and_summarize_opt(
        atoms,
        preset=preset,
        calc_defaults=calc_defaults,
        calc_swaps=calc_kwargs,
        opt_defaults=opt_defaults,
        opt_params=opt_params,
        additional_fields={"name": "VASP ASE Relax"},
        copy_files=copy_files,
    )


@job
def non_scf_job(
    atoms: Atoms,
    prev_dir: SourceDirectory,
    preset: str | None = "BulkSet",
    nbands_factor: float = 1.2,
    kpts_mode: Literal["uniform", "line"] = "uniform",
    uniform_kppvol: float = 100,
    line_kpt_density: float = 20,
    calculate_optics: bool = False,
    **calc_kwargs,
) -> VaspSchema:
    """
    Carry out a non-self-consistent field (NSCF) calculation.

    Parameters
    ----------
    atoms
        Atoms object.
    prev_dir
        Directory of the prior job. Must contain a CHGCAR and vasprun.xml file.
    preset
        Preset to use from `quacc.calculators.vasp.presets`.
    nbands_factor
        A multiplicative factor used to adjust NBANDS when vasprun.xml(.gz) exists in
        prev_dir
    kpts_mode
        Type of k-points mode. Options are "uniform" or "line".
    uniform_kppvol
        The k-point per volume density for the uniform k-point mode.
    line_kpt_density
        The k-point density for the line k-point mode.
    calculate_optics
        Whether to calculate optical properties.
    **calc_kwargs
        Custom kwargs for the Vasp calculator. Set a value to
        `None` to remove a pre-existing key entirely. For a list of available
        keys, refer to [quacc.calculators.vasp.vasp.Vasp][].

    Returns
    -------
    VaspSchema
        Dictionary of results from [quacc.schemas.vasp.VaspSummarize.run][].
        See the type-hint for the data structure.
    """

    vasprun_path = zpath(str(Path(prev_dir, "vasprun.xml")))
    vasprun = Vasprun(vasprun_path)

    prior_nbands = vasprun.parameters["NBANDS"]
    calc_defaults: dict[str, Any] = {
        "icharg": 11,
        "kspacing": None,
        "lcharg": False,
        "lorbit": 11,
        "lwave": False,
        "nsw": 0,
        "nbands": int(np.ceil(prior_nbands * nbands_factor)),
    }
    if kpts_mode == "uniform":
        calc_defaults |= {
            "ismear": -5,
            "isym": 2,
            "pmg_kpts": {"kppvol": uniform_kppvol},
            "nedos": 6001,
        }
    elif kpts_mode == "line":
        is_metal = vasprun.get_band_structure().is_metal()
        calc_defaults |= {
            "ismear": 1 if is_metal else 0,
            "isym": 0,
            "pmg_kpts": {"line_density": line_kpt_density},
            "sigma": 0.2 if is_metal else 0.01,
        }
    else:
        raise ValueError("Supported kpoint modes are 'uniform' and 'line' at present")

    if calculate_optics:
        calc_defaults |= {"cshift": 1e-5, "loptics": True, "lreal": False}

    return run_and_summarize(
        atoms,
        preset=preset,
        calc_defaults=calc_defaults,
        calc_swaps=calc_kwargs,
        additional_fields={"name": "VASP Non-SCF"},
        copy_files={prev_dir: ["CHGCAR*", "WAVECAR*"]},
    )


@job
def freq_job(
    atoms: Atoms,
    preset: str | None = "BulkSet",
    energy: float = 0.0,
    temperature: float = 298.15,
    pressure: float = 1.0,
    thermo_method: Literal["harmonic", "ideal_gas"] = "harmonic",
    vib_kwargs: VibKwargs | None = None,
    copy_files: SourceDirectory | dict[SourceDirectory, Filenames] | None = None,
    **calc_kwargs,
) -> VibThermoSchema:
    """
    Run a frequency job and calculate thermochemistry.

    Parameters
    ----------
    atoms
        Atoms object
    preset
        Preset to use from `quacc.calculators.vasp.presets`.
    energy
        Potential energy in eV. If 0, then the output is just the correction.
    temperature
        Temperature in Kelvins.
    pressure
        Pressure in bar.
    thermo_method
        Method to use for thermochemistry. Options are "harmonic" or "ideal_gas".
    vib_kwargs
        Dictionary of kwargs for the [ase.vibrations.Vibrations][] class.
    copy_files
        Files to copy (and decompress) from source to the runtime directory.
    **calc_kwargs
        Custom kwargs for the Vasp calculator. Set a value to
        `None` to remove a pre-existing key entirely. For a list of available
        keys, refer to [quacc.calculators.vasp.vasp.Vasp][].

    Returns
    -------
    VibThermoSchema
        Dictionary of results
    """
    calc_defaults = {"ediff": 1e-7, "isym": 0, "lcharg": False, "lwave": True, "nsw": 0}
    vib_kwargs = vib_kwargs or {}

    return run_and_summarize_vib_and_thermo(
        atoms,
        energy=energy,
        temperature=temperature,
        pressure=pressure,
        thermo_method=thermo_method,
        preset=preset,
        calc_defaults=calc_defaults,
        calc_swaps=calc_kwargs,
        vib_kwargs=vib_kwargs,
        copy_files=copy_files,
        additional_fields={"name": "VASP Frequency and Thermo"},
    )

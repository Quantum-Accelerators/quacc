"""Core recipes for the tblite code."""
from __future__ import annotations

from typing import TYPE_CHECKING

from ase.optimize import FIRE
from monty.dev import requires

from quacc import job
from quacc.builders.thermo import build_ideal_gas
from quacc.runners.calc import run_ase_calc, run_ase_opt, run_ase_phonons, run_ase_vib
from quacc.schemas.ase import (
    summarize_opt_run,
    summarize_phonon_run,
    summarize_run,
    summarize_vib_and_thermo,
)
from quacc.utils.dicts import merge_dicts

try:
    from tblite.ase import TBLite
except ImportError:
    TBLite = None

if TYPE_CHECKING:
    from typing import Any, Literal

    from ase import Atoms

    from quacc.runners.calc import PhononKwargs, PhononReadKwargs, VibKwargs
    from quacc.schemas._aliases.ase import (
        OptSchema,
        PhononSchema,
        RunSchema,
        VibThermoSchema,
    )


@job
@requires(TBLite, "tblite must be installed. Refer to the quacc documentation.")
def static_job(
    atoms: Atoms,
    method: Literal["GFN1-xTB", "GFN2-xTB", "IPEA1-xTB"] = "GFN2-xTB",
    calc_swaps: dict[str, Any] | None = None,
    copy_files: list[str] | None = None,
) -> RunSchema:
    """
    Carry out a single-point calculation.

    Parameters
    ----------
    atoms
        Atoms object
    method
        GFN1-xTB, GFN2-xTB, and IPEA1-xTB.
    calc_swaps
        Dictionary of custom kwargs for the EMT calculator. Set a value to
        `None` to remove a pre-existing key entirely. For a list of available
        keys, refer to the `tblite.ase.TBLite` calculator.

        !!! Info "Calculator defaults"

            ```python
            {"method": method}
            ```
    copy_files
        Files to copy to the runtime directory.

    Returns
    -------
    RunSchema
        Dictionary of results from [quacc.schemas.ase.summarize_run][]
    """

    defaults = {"method": method}
    flags = merge_dicts(defaults, calc_swaps)
    atoms.calc = TBLite(**flags)

    final_atoms = run_ase_calc(atoms, copy_files=copy_files)
    return summarize_run(
        final_atoms,
        input_atoms=atoms,
        additional_fields={"name": "TBLite Static"},
    )


@job
@requires(TBLite, "tblite must be installed. Refer to the quacc documentation.")
def relax_job(
    atoms: Atoms,
    method: Literal["GFN1-xTB", "GFN2-xTB", "IPEA1-xTB"] = "GFN2-xTB",
    relax_cell: bool = False,
    calc_swaps: dict[str, Any] | None = None,
    opt_swaps: dict[str, Any] | None = None,
    copy_files: list[str] | None = None,
) -> OptSchema:
    """
    Relax a structure.

    !!! Info "Calculator defaults"

        ```python
        {"method": method}
        ```

    !!! Info "Optimizer defaults"

        ```python
        {"fmax": 0.01, "max_steps": 1000, "optimizer": FIRE}
        ```

    Parameters
    ----------
    atoms
        Atoms object
    method
        GFN0-xTB, GFN1-xTB, GFN2-xTB.
    relax_cell
        Whether to relax the cell.
    calc_swaps
        Dictionary of custom kwargs for the tblite calculator. Set a value to
        `None` to remove a pre-existing key entirely. For a list of available
        keys, refer to the `tblite.ase.TBLite` calculator.
    opt_swaps
        Dictionary of custom kwargs for the optimization process. Set a value
        to `None` to remove a pre-existing key entirely. For a list of available
        keys, refer to [quacc.runners.calc.run_ase_opt][].
    copy_files
        Files to copy to the runtime directory.

    Returns
    -------
    OptSchema
        Dictionary of results from [quacc.schemas.ase.summarize_opt_run][]
    """

    defaults = {"method": method}
    flags = merge_dicts(defaults, calc_swaps)
    atoms.calc = TBLite(**flags)

    opt_defaults = {"fmax": 0.01, "max_steps": 1000, "optimizer": FIRE}
    opt_flags = merge_dicts(opt_defaults, opt_swaps)

    dyn = run_ase_opt(atoms, relax_cell=relax_cell, copy_files=copy_files, **opt_flags)

    return summarize_opt_run(dyn, additional_fields={"name": "TBLite Relax"})


@job
@requires(TBLite, "tblite must be installed. Refer to the quacc documentation.")
def freq_job(
    atoms: Atoms,
    method: Literal["GFN1-xTB", "GFN2-xTB", "IPEA1-xTB"] = "GFN2-xTB",
    energy: float = 0.0,
    temperature: float = 298.15,
    pressure: float = 1.0,
    calc_swaps: dict[str, Any] | None = None,
    vib_kwargs: VibKwargs | None = None,
    copy_files: list[str] | None = None,
) -> VibThermoSchema:
    """
    Run a frequency job and calculate thermochemistry.

    !!! Info "Calculator defaults"

        ```python
        {"method": method}
        ```

    Parameters
    ----------
    atoms
        Atoms object
    method
        GFN0-xTB, GFN1-xTB, GFN2-xTB, GFN-FF.
    energy
        Potential energy in eV. If 0, then the output is just the correction.
    temperature
        Temperature in Kelvins.
    pressure
        Pressure in bar.
    calc_swaps
        Dictionary of custom kwargs for the tblite calculator. Set a value to
        `None` to remove a pre-existing key entirely. For a list of available
        keys, refer to the `tblite.ase.TBLite` calculator.
    vib_kwargs
        Dictionary of custom kwargs for the vibration analysis. Refer to
        [quacc.runners.calc.run_ase_vib][].
    copy_files
        Files to copy to the runtime directory.

    Returns
    -------
    VibThermoSchema
        Dictionary of results from [quacc.schemas.ase.summarize_vib_and_thermo][]
    """
    vib_kwargs = vib_kwargs or {}

    defaults = {"method": method}
    flags = merge_dicts(defaults, calc_swaps)
    atoms.calc = TBLite(**flags)

    vibrations = run_ase_vib(atoms, vib_kwargs=vib_kwargs, copy_files=copy_files)
    igt = build_ideal_gas(atoms, vibrations.get_frequencies(), energy=energy)

    return summarize_vib_and_thermo(
        vibrations,
        igt,
        temperature=temperature,
        pressure=pressure,
        additional_fields={"name": "TBLite Frequency and Thermo"},
    )


@job
@requires(TBLite, "tblite must be installed. Refer to the quacc documentation.")
def phonon_job(
    atoms: Atoms,
    method: Literal["GFN1-xTB", "GFN2-xTB", "IPEA1-xTB"] = "GFN2-xTB",
    calc_swaps: dict[str, Any] | None = None,
    phonon_kwargs: PhononKwargs | PhononReadKwargs | None = None,
    copy_files: list[str] | None = None,
) -> PhononSchema:
    """
    Relax a structure.

    !!! Info "Calculator defaults"

        ```python
        {"method": method}
        ```

    Parameters
    ----------
    atoms
        Atoms object
    method
        GFN0-xTB, GFN1-xTB, GFN2-xTB.
    relax_cell
        Whether to relax the cell.
    calc_swaps
        Dictionary of custom kwargs for the tblite calculator. Set a value to
        `None` to remove a pre-existing key entirely. For a list of available
        keys, refer to the `tblite.ase.TBLite` calculator.
    phonon_kwargs
        Dictionary of custom kwargs for [quacc.runners.calc.run_ase_vib][]
    copy_files
        Files to copy to the runtime directory.

    Returns
    -------
    PhononSchema
        Dictionary of results from [quacc.schemas.ase.summarize_phonon_run][]
    """
    phonon_kwargs = phonon_kwargs or {}

    defaults = {"method": method}
    flags = merge_dicts(defaults, calc_swaps)
    atoms.calc = TBLite(**flags)

    dyn = run_ase_phonons(atoms, copy_files=copy_files, **phonon_kwargs)

    return summarize_phonon_run(dyn, additional_fields={"name": "TBLite Phonons"})

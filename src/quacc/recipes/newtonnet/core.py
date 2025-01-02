"""Core recipes for the NewtonNet code."""

from __future__ import annotations

from importlib.util import find_spec
from typing import TYPE_CHECKING

from ase.vibrations.data import VibrationsData
from monty.dev import requires

from quacc import get_settings, job
from quacc.runners.ase import Runner
from quacc.schemas.ase import Summarize, VibSummarize
from quacc.utils.dicts import recursive_dict_merge

has_sella = bool(find_spec("sella"))
has_newtonnet = bool(find_spec("newtonnet"))

if has_sella:
    from sella import Sella
if has_newtonnet:
    from newtonnet.utils.ase_interface import MLAseCalculator as NewtonNet

if TYPE_CHECKING:
    from typing import Any

    from ase.atoms import Atoms

    from quacc.types import (
        Filenames,
        OptParams,
        OptSchema,
        RunSchema,
        SourceDirectory,
        VibThermoSchema,
    )


@job
@requires(
    has_newtonnet, "NewtonNet must be installed. Refer to the quacc documentation."
)
def static_job(
    atoms: Atoms,
    copy_files: SourceDirectory | dict[SourceDirectory, Filenames] | None = None,
    additional_fields: dict[str, Any] | None = None,
    **calc_kwargs,
) -> RunSchema:
    """
    Carry out a single-point calculation.

    Parameters
    ----------
    atoms
        Atoms object
    copy_files
        Files to copy (and decompress) from source to the runtime directory.
    additional_fields
        Additional fields to add to the results dictionary.
    **calc_kwargs
        Custom kwargs for the NewtonNet calculator. Set a value to
        `quacc.Remove` to remove a pre-existing key entirely. For a list of available
        keys, refer to the `newtonnet.utils.ase_interface.MLAseCalculator` calculator.

    Returns
    -------
    RunSchema
        Dictionary of results, specified in [quacc.schemas.ase.Summarize.run][].
        See the type-hint for the data structure.
    """
    settings = get_settings()
    calc_defaults = {
        "model_path": settings.NEWTONNET_MODEL_PATH,
        "settings_path": settings.NEWTONNET_CONFIG_PATH,
    }
    calc_flags = recursive_dict_merge(calc_defaults, calc_kwargs)

    calc = NewtonNet(**calc_flags)
    final_atoms = Runner(atoms, calc, copy_files=copy_files).run_calc()

    return Summarize(
        additional_fields={"name": "NewtonNet Static"} | (additional_fields or {})
    ).run(final_atoms, atoms)


@job
@requires(
    has_newtonnet, "NewtonNet must be installed. Refer to the quacc documentation."
)
def relax_job(
    atoms: Atoms,
    opt_params: OptParams | None = None,
    copy_files: SourceDirectory | dict[SourceDirectory, Filenames] | None = None,
    additional_fields: dict[str, Any] | None = None,
    **calc_kwargs,
) -> OptSchema:
    """
    Relax a structure.

    Parameters
    ----------
    atoms
        Atoms object
    opt_params
        Dictionary of custom kwargs for the optimization process. For a list
        of available keys, refer to [quacc.runners.ase.Runner.run_opt][].
    copy_files
        Files to copy (and decompress) from source to the runtime directory.
    additional_fields
        Additional fields to add to the results dictionary.
    **calc_kwargs
        Dictionary of custom kwargs for the NewtonNet calculator. Set a value to
        `quacc.Remove` to remove a pre-existing key entirely. For a list of available
        keys, refer to the `newtonnet.utils.ase_interface.MLAseCalculator` calculator.

    Returns
    -------
    OptSchema
        Dictionary of results, specified in [quacc.schemas.ase.Summarize.opt][].
        See the type-hint for the data structure.
    """
    settings = get_settings()
    calc_defaults = {
        "model_path": settings.NEWTONNET_MODEL_PATH,
        "settings_path": settings.NEWTONNET_CONFIG_PATH,
    }
    opt_defaults = {"optimizer": Sella} if has_sella else {}

    calc_flags = recursive_dict_merge(calc_defaults, calc_kwargs)
    opt_flags = recursive_dict_merge(opt_defaults, opt_params)

    calc = NewtonNet(**calc_flags)
    dyn = Runner(atoms, calc, copy_files=copy_files).run_opt(**opt_flags)

    return _add_stdev_and_hess(
        Summarize(
            additional_fields={"name": "NewtonNet Relax"} | (additional_fields or {})
        ).opt(dyn)
    )


@job
@requires(
    has_newtonnet, "NewtonNet must be installed. Refer to the quacc documentation."
)
def freq_job(
    atoms: Atoms,
    temperature: float = 298.15,
    pressure: float = 1.0,
    copy_files: SourceDirectory | dict[SourceDirectory, Filenames] | None = None,
    additional_fields: dict[str, Any] | None = None,
    **calc_kwargs,
) -> VibThermoSchema:
    """
    Perform a frequency calculation using the given atoms object.

    Parameters
    ----------
    atoms
        The atoms object representing the system.
    temperature
        The temperature for the thermodynamic analysis.
    pressure
        The pressure for the thermodynamic analysis.
    copy_files
        Files to copy (and decompress) from source to the runtime directory.
    additional_fields
        Additional fields to add to the results dictionary.
    **calc_kwargs
        Custom kwargs for the NewtonNet calculator. Set a value to
        `quacc.Remove` to remove a pre-existing key entirely. For a list of available
        keys, refer to the `newtonnet.utils.ase_interface.MLAseCalculator` calculator.

    Returns
    -------
    VibThermoSchema
        Dictionary of results. See the type-hint for the data structure.
    """
    settings = get_settings()
    calc_defaults = {
        "model_path": settings.NEWTONNET_MODEL_PATH,
        "settings_path": settings.NEWTONNET_CONFIG_PATH,
        "hess_method": "autograd",
    }
    calc_flags = recursive_dict_merge(calc_defaults, calc_kwargs)

    calc = NewtonNet(**calc_flags)
    final_atoms = Runner(atoms, calc, copy_files=copy_files).run_calc()

    summary = Summarize(
        additional_fields={"name": "NewtonNet Frequency"} | (additional_fields or {})
    ).run(final_atoms, atoms)

    vib = VibrationsData(final_atoms, summary["results"]["hessian"])
    return VibSummarize(
        vib,
        directory=summary["dir_name"],
        additional_fields={"name": "ASE Vibrations and Thermo Analysis"},
    ).vib_and_thermo(
        "ideal_gas",
        energy=summary["results"]["energy"],
        temperature=temperature,
        pressure=pressure,
    )


def _add_stdev_and_hess(summary: dict[str, Any], **calc_kwargs) -> dict[str, Any]:
    """
    Calculate and add standard deviation values and Hessians to the summary.

    This function takes a summary dictionary containing information about a
    molecular trajectory and calculates the standard deviation of various
    properties using the NewtonNet machine learning calculator. It adds the
    calculated standard deviation values and Hessians to each configuration in
    the trajectory.

    Parameters
    ----------
    summary
        A dictionary containing information about the molecular trajectory.
    **calc_kwargs
        Custom kwargs for the NewtonNet calculator. Set a value to
        `quacc.Remove` to remove a pre-existing key entirely. For a list of available
        keys, refer to the `newtonnet.utils.ase_interface.MLAseCalculator` calculator.

    Returns
    -------
    dict[str, Any]
        The modified summary dictionary with added standard deviation and
        Hessian values.
    """
    settings = get_settings()
    calc_defaults = {
        "model_path": settings.NEWTONNET_MODEL_PATH,
        "settings_path": settings.NEWTONNET_CONFIG_PATH,
    }
    calc_flags = recursive_dict_merge(calc_defaults, calc_kwargs)
    for i, atoms in enumerate(summary["trajectory"]):
        calc = NewtonNet(**calc_flags)
        results = Runner(atoms, calc).run_calc().calc.results
        summary["trajectory_results"][i]["hessian"] = results["hessian"]
        summary["trajectory_results"][i]["energy_std"] = results["energy_disagreement"]
        summary["trajectory_results"][i]["forces_std"] = results["forces_disagreement"]
        summary["trajectory_results"][i]["hessian_std"] = results[
            "hessian_disagreement"
        ]

    return summary

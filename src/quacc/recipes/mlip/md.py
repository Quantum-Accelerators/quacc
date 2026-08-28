"""Molecular dynamics recipes for universal machine-learned interatomic potentials."""

from __future__ import annotations

from typing import TYPE_CHECKING

from ase.md.npt import NPT
from ase.md.verlet import VelocityVerlet
from ase.units import bar, fs

from quacc import Remove, job
from quacc.recipes.mlip._base import pick_calculator
from quacc.runners.ase import Runner
from quacc.schemas.ase import Summarize
from quacc.utils.dicts import recursive_dict_merge
from quacc.utils.md import resolve_md_ensemble, upper_triangular_cell

if TYPE_CHECKING:
    from typing import Any, Literal

    from ase.atoms import Atoms
    from ase.md.md import MolecularDynamics

    from quacc.types import DynSchema, MDEnsemble, MDParams


@job
def md_job(
    atoms: Atoms,
    library: Literal["fairchem", "matcalc", "rootstock"],
    dynamics: MolecularDynamics | MDEnsemble = VelocityVerlet,
    steps: int = 10000,
    timestep_fs: float = 0.5,
    temperature_K: float | None = None,
    pressure_bar: float | None = None,
    md_params: MDParams | None = None,
    additional_fields: dict[str, Any] | None = None,
    **calc_kwargs: Any,
) -> DynSchema:
    """
    Carry out a Molecular Dynamics calculation.

    Parameters
    ----------
    atoms
        Atoms object
    library
        MLIP library to use:
        - `fairchem` passes `**calc_kwargs` to `FAIRChemCalculator.from_model_checkpoint()`
        - `matcalc` passes `**calc_kwargs` to `matcalc.load_fp()`
        - `rootstock` passes `**calc_kwargs` to `rootstock.RootstockCalculator()`
    dynamics
        Either an ASE `MolecularDynamics` class (from `ase.md.md.MolecularDynamics`)
        or an ensemble name. Supported ensemble names, their ASE dynamics classes,
        and their default parameters (any of which can be overridden via
        `md_params["dynamics_kwargs"]`; `taut` defaults to 100x the timestep and
        `taup` to 1000x the timestep):

        - `"nve"`: `VelocityVerlet` (`temperature_K` and `pressure_bar` are ignored)
        - `"nvt"`: `NoseHooverChainNVT` (`tdamp=taut`)
        - `"nvt_berendsen"`: `NVTBerendsen` (`taut`)
        - `"nvt_langevin"`: `Langevin` (`friction=0.001`)
        - `"nvt_andersen"`: `Andersen` (`andersen_prob=0.01`)
        - `"nvt_bussi"`: `Bussi` (`taut`)
        - `"npt"`: `NPT` (`ttime=25 fs`, `pfactor=75**2 fs`)
        - `"npt_berendsen"`: `NPTBerendsen` (`taut`, `taup`)
        - `"npt_inhomogeneous"`: `Inhomogeneous_NPTBerendsen` (`taut`, `taup`)
        - `"npt_mtk"`: `MTKNPT` (`tdamp=taut`, `pdamp=taup`)
        - `"npt_isotropic_mtk"`: `IsotropicMTKNPT` (`tdamp=taut`, `pdamp=taup`)

        For the NPT-family ensembles, `pressure_bar` is routed to the appropriate
        keyword argument of the dynamics class (`externalstress` or `pressure_au`)
        and defaults to 0 bar when not specified. If the `NPT` class is used, the
        cell is transformed to the upper-triangular form it requires.
    steps
        Number of MD steps to run.
    timestep_fs
        Time step in fs.
    temperature_K
        Temperature in K, if applicable for the given ensemble.
    pressure_bar
        Pressure in bar, if applicable for the given ensemble. When `dynamics` is
        an ASE class, this is passed as `pressure_au` and is only compatible with
        classes that accept that keyword argument (e.g. the Berendsen family);
        prefer an ensemble name for NPT runs.
    md_params
        Dictionary of custom kwargs for the MD run. For a list of available
        keys, refer to [quacc.runners.ase.Runner.run_md][].
    additional_fields
        Additional fields to add to the results dictionary.
    **calc_kwargs
        Custom kwargs for the underlying MLIP library.

    Returns
    -------
    DynSchema
        Dictionary of results, specified in [quacc.schemas.ase.Summarize.md][].
        See the type-hint for the data structure.
    """
    if isinstance(dynamics, str):
        dynamics, dynamics_defaults = resolve_md_ensemble(
            dynamics, timestep_fs * fs, temperature_K, pressure_bar
        )
    else:
        dynamics_defaults = {
            "temperature_K": temperature_K if temperature_K is not None else Remove,
            "pressure_au": pressure_bar * bar if pressure_bar is not None else Remove,
        }
    dynamics_defaults["timestep"] = timestep_fs * fs

    if dynamics is NPT:
        atoms = upper_triangular_cell(atoms)

    md_defaults = {"steps": steps, "dynamics_kwargs": dynamics_defaults}
    md_params = recursive_dict_merge(md_defaults, md_params)

    calc = pick_calculator(library, **calc_kwargs)
    dyn = Runner(atoms, calc).run_md(dynamics, **md_params)

    return Summarize(
        additional_fields={"name": f"{library} MD"} | (additional_fields or {})
    ).md(dyn)

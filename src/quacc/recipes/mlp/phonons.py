"""Phonon recipes for MLPs."""
from __future__ import annotations

from typing import TYPE_CHECKING

from quacc import flow
from quacc.recipes.common.phonons import phonon_flow as common_phonon_flow
from quacc.recipes.mlp.core import relax_job, static_job
from quacc.utils.dicts import recursive_dict_merge
from quacc.wflow_tools.customizers import customize_funcs

if TYPE_CHECKING:
    from typing import Any, Callable, Literal

    from ase.atoms import Atoms

    from quacc.schemas._aliases.phonons import PhononSchema


@flow
def phonon_flow(
    atoms: Atoms,
    method: Literal["mace", "m3gnet", "chgnet"],
    symprec: float = 1e-4,
    min_length: float | None = 15.0,
    displacement: float = 0.01,
    t_step: float = 10,
    t_min: float = 0,
    t_max: float = 1000,
    run_relax: bool = True,
    job_params: dict[str, dict[str, Any]] | None = None,
    job_decorators: dict[str, Callable | None] | None = None,
) -> PhononSchema:
    """
    Carry out a phonon workflow, consisting of:

    1. Optional relaxation.
        - name: "relax_job"
        - job: [quacc.recipes.mlp.core.relax_job][]

    2. Generation of supercells.

    3. Static calculations on supercells
        - name: "static_job"
        - job: [quacc.recipes.mlp.core.static_job][]

    4. Calculation of thermodynamic properties.

    Parameters
    ----------
    atoms
        Atoms object
    method
        Universal ML interatomic potential method to use
    symprec
        Precision for symmetry detection.
    min_length
        Minimum length of each lattice dimension (A).
    displacement
        Atomic displacement (A).
    t_step
        Temperature step (K).
    t_min
        Min temperature (K).
    t_max
        Max temperature (K).
    job_params
        Custom parameters to pass to each Job in the Flow. This is a dictinoary where
        the keys are the names of the jobs and the values are dictionaries of parameters.
    job_decorators
        Custom decorators to apply to each Job in the Flow. This is a dictionary where
        the keys are the names of the jobs and the values are decorators.

    Returns
    -------
    PhononSchema
        Dictionary of results from [quacc.schemas.phonons.summarize_phonopy][].
        See the type-hint for the data structure.
    """
    calc_defaults = {
        "relax_job": {"method": method, "opt_params": {"fmax": 1e-3}},
        "static_job": {"method": method},
    }
    job_params = recursive_dict_merge(calc_defaults, job_params)

    relax_job_, static_job_ = customize_funcs(
        ["relax_job", "static_job"],
        [relax_job, static_job],
        parameters=job_params,
        decorators=job_decorators,
    )

    return common_phonon_flow(
        atoms,
        static_job_,
        relax_job=relax_job_ if run_relax else None,
        symprec=symprec,
        min_length=min_length,
        displacement=displacement,
        t_step=t_step,
        t_min=t_min,
        t_max=t_max,
        additional_fields={"name": f"{method} Phonons"},
    )

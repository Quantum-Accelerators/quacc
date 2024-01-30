"""Phonon recipes for EMT."""
from __future__ import annotations

from typing import TYPE_CHECKING

from quacc import flow
from quacc.recipes.common.phonons import phonon_flow as common_phonon_flow
from quacc.recipes.emt.core import static_job
from quacc.wflow_tools.customizers import customize_funcs

if TYPE_CHECKING:
    from typing import Any, Callable

    from ase.atoms import Atoms

    from quacc.schemas._aliases.phonons import PhononSchema


@flow
def phonon_flow(
    atoms: Atoms,
    symprec: float = 1e-4,
    min_length: float | None = 20.0,
    atom_disp: float = 0.01,
    t_step: float = 10,
    t_min: float = 0,
    t_max: float = 1000,
    job_params: dict[str, dict[str, Any]] | None = None,
    job_decorators: dict[str, Callable | None] | None = None,
) -> PhononSchema:
    """
    Carry out a phonon workflow, consisting of:

    1. Generation of supercells.

    2. Static calculations on supercells
        - name: "static_job"
        - job: [quacc.recipes.emt.core.static_job][]

    3. Calculation of thermodynamic properties.

    Parameters
    ----------
    atoms
        Atoms object
    symprec
        Precision for symmetry detection.
    min_length
        Minimum length of each lattice dimension (A).
    atom_disp
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
        See the return type-hint for the data structure.
    """
    static_job_ = customize_funcs(
        "static_job", static_job, parameters=job_params, decorators=job_decorators
    )

    return common_phonon_flow(
        atoms,
        static_job_,
        symprec=symprec,
        min_length=min_length,
        atom_disp=atom_disp,
        t_step=t_step,
        t_min=t_min,
        t_max=t_max,
        additional_fields={"name": "EMT Phonons"},
    )

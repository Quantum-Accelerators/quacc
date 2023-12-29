"""Phonon recipes for EMT"""
from __future__ import annotations

from typing import TYPE_CHECKING

from quacc import flow
from quacc.recipes.common.phonons import phonon_flow as phonon_flow_
from quacc.recipes.emt.core import static_job
from quacc.wflow_tools.customizers import customize_funcs

if TYPE_CHECKING:
    from typing import Any, Callable

    from ase.atoms import Atoms
    from numpy.typing import ArrayLike

    from quacc.schemas._aliases.phonons import PhononSchema


@flow
def phonon_flow(
    atoms: Atoms,
    supercell_matrix: ArrayLike = ((2, 0, 0), (0, 2, 0), (0, 0, 2)),
    atom_disp: float = 0.01,
    symprec: float = 1e-5,
    t_step: float = 10,
    t_min: float = 0,
    t_max: float = 1000,
    decorators: dict[str, Callable | None] | None = None,
    parameters: dict[str, Any] | None = None,
) -> PhononSchema:
    """
    Carry out a phonon calculation.

    Parameters
    ----------
    atoms
        Atoms object
    supercell_matrix
        Supercell matrix to use. Defaults to 2x2x2 supercell.
    atom_disp
        Atomic displacement (A).
    symprec
        Precision for symmetry detection.
    t_step
        Temperature step (K).
    t_min
        Min temperature (K).
    t_max
        Max temperature (K).
    decorators
        Custom decorators to apply to each Job in the Flow.
        Refer to [quacc.wflow_tools.customizers.customize_funcs][] for details.
    parameters
        Custom parameters to pass to each Job in the Flow.
        Refer to [quacc.wflow_tools.customizers.customize_funcs][] for details.

    Returns
    -------
    PhononSchema
        Dictionary of results from [quacc.schemas.phonons.summarize_phonopy][]
    """
    static_job_ = customize_funcs(
        {"static_job": static_job}, decorators=decorators, parameters=parameters
    )[0]

    return phonon_flow_(
        atoms,
        static_job_,
        supercell_matrix=supercell_matrix,
        atom_disp=atom_disp,
        t_step=t_step,
        t_min=t_min,
        t_max=t_max,
        phonopy_kwargs={"symprec": symprec},
        additional_fields={"name": "EMT Phonons"},
    )

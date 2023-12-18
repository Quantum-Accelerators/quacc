"""Phonon recipes for TBLite"""
from __future__ import annotations

from typing import TYPE_CHECKING

from quacc import flow
from quacc.recipes.common.phonons import phonon_flow as phonon_flow_
from quacc.recipes.tblite.core import static_job as static_job_

if TYPE_CHECKING:
    from ase.atoms import Atoms
    from numpy.typing import ArrayLike

    from quacc import Job
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
    static_job: Job | None = static_job_,
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
    static_job
        Static job, which defaults to [quacc.recipes.tblite.core.static_job][]

    Returns
    -------
    PhononSchema
        Dictionary of results from [quacc.schemas.phonons.summarize_phonopy][]
    """

    return phonon_flow_(
        atoms,
        static_job,
        supercell_matrix=supercell_matrix,
        atom_disp=atom_disp,
        t_step=t_step,
        t_min=t_min,
        t_max=t_max,
        phonopy_kwargs={"symprec": symprec},
        additional_fields={"name": "TBLite Phonons"},
    )

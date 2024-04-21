"""Phonon recipes for MLPs."""

from __future__ import annotations

from importlib.util import find_spec
from typing import TYPE_CHECKING

from monty.dev import requires

from quacc import flow
from quacc.recipes.common.phonons import phonon_flow as common_phonon_flow
from quacc.recipes.mlp.core import relax_job, static_job
from quacc.utils.dicts import recursive_dict_merge
from quacc.wflow_tools.customizers import customize_funcs

has_deps = find_spec("phonopy") is not None and find_spec("seekpath") is not None

if TYPE_CHECKING:
    from typing import Any, Callable, Literal

    from ase.atoms import Atoms

    from quacc.schemas._aliases.phonons import PhononSchema


@flow
@requires(
    has_deps,
    message="Phonopy and seekpath must be installed. Run `pip install quacc[phonons]`",
)
def phonon_flow(
    atoms: Atoms,
    method: Literal["mace", "m3gnet", "chgnet"],
    fixed_atoms: list[int] | None = None,
    min_lengths: float | tuple[float, float, float] | None = 20.0,
    phonopy_kwargs: dict[str, Any] | None = None,
    force_constants_kwargs: dict[str, Any] | None = None,
    mesh_kwargs: dict[str, Any] | None = None,
    total_dos_kwargs: dict[str, Any] | None = None,
    thermal_properties_kwargs: dict[str, Any] | None = None,
    auto_band_structure_kwargs: dict[str, Any] | None = None,
    generate_displacements_kwargs: dict[str, Any] | None = None,
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
    min_lengths
        Minimum length of each lattice dimension (A).
    supercell_matrix
        The supercell matrix to use. If specified, it will override any
        value specified by `min_lengths`.
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
        fixed_atoms=fixed_atoms,
        min_lengths=min_lengths,
        phonopy_kwargs=phonopy_kwargs,
        force_constants_kwargs=force_constants_kwargs,
        mesh_kwargs=mesh_kwargs,
        total_dos_kwargs=total_dos_kwargs,
        thermal_properties_kwargs=thermal_properties_kwargs,
        auto_band_structure_kwargs=auto_band_structure_kwargs,
        generate_displacements_kwargs=generate_displacements_kwargs,
        additional_fields={"name": f"{method} Phonons"},
    )

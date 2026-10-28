"""Common equation of state workflows."""

from __future__ import annotations

from typing import TYPE_CHECKING

from emmet.core.eos import EOSDoc

from quacc import flow, job, subflow
from quacc.atoms.deformation import make_volume_strains_from_bulk

if TYPE_CHECKING:
    from collections.abc import Callable
    from typing import Any, TypeVar

    from ase.atoms import Atoms

    from quacc.types import EOSSchema, RunSchema

    Schema = TypeVar("Schema", bound=RunSchema)


@flow
def eos_flow(
    atoms: Atoms,
    relax_job: Callable[..., Schema],
    static_job: Callable[..., Schema],
    pre_relax: bool = True,
    run_static: bool = False,
    strain_kwargs: dict[str, Any] | None = None,
    include_undeformed: bool = True,
) -> EOSSchema[Schema]:
    """
    Common workflow for fitting an equation of state.

    Parameters
    ----------
    atoms
        Atoms object
    relax_job
        The relaxation function. It is called with `relax_cell=True` on the bulk
        structure and with `relax_cell=False` on each volume-strained structure.
    static_job
        The static function.
    pre_relax
        Whether to run a cell relaxation on the bulk structure before straining it
        (true) or run a static calculation (false).
    run_static
        Whether to run static calculations after the relaxations on the undeformed
        and volume-strained structures. Recommended when the final energy of a
        relaxation is less accurate than a static calculation (e.g. VASP), and
        unnecessary when it is not (e.g. EMT or MLIPs).
    strain_kwargs
        Additional keyword arguments to pass to
        [quacc.atoms.deformation.make_volume_strains_from_bulk][]
    include_undeformed
        Whether to include the energy of the undeformed result as a data point in
        the fit. See [quacc.recipes.common.eos.eos_subflow][].

    Returns
    -------
    EOSSchema
        See the return type-hint for the data structure.
    """
    if pre_relax:
        undeformed_result = relax_job(atoms, relax_cell=True)
        if run_static:
            undeformed_result = static_job(undeformed_result["atoms"])
    else:
        undeformed_result = static_job(atoms)

    return eos_subflow(
        undeformed_result,
        relax_job=relax_job,
        static_job=static_job if run_static else None,
        strain_kwargs=strain_kwargs,
        include_undeformed=include_undeformed,
    )


@subflow
def eos_subflow(
    undeformed_result: RunSchema,
    relax_job: Callable[..., Schema] | None = None,
    static_job: Callable[..., Schema] | None = None,
    strain_kwargs: dict[str, Any] | None = None,
    include_undeformed: bool = True,
) -> EOSSchema[Schema]:
    """
    Workflow consisting of:

    1. Volume-strained structures generation

    2. Volume-strained structures relaxations at fixed cell (optional)

    3. Volume-strained structures statics (optional)

    4. Equation of state fit

    This workflow takes the result of a prior calculation (typically a relaxation
    and/or static calculation on the bulk structure) and can be used with any
    calculator by supplying the corresponding jobs. Calculator-specific settings
    should be applied to the jobs before they are passed in, e.g. with
    [quacc.wflow_tools.customizers.update_parameters][].

    At least one of `relax_job` and `static_job` must be supplied. The fit uses the
    static energies if `static_job` is supplied and the relaxation energies
    otherwise. Supplying only `static_job` keeps the fractional coordinates fixed,
    which is only exact when symmetry fixes all atomic positions (e.g. fcc metals).

    Parameters
    ----------
    undeformed_result
        Result of a static or optimization calculation on the bulk structure.
    relax_job
        The relaxation function to run on each volume-strained structure. It is
        called with `relax_cell=False`, so only the atomic positions relax.
    static_job
        The static function to run on each volume-strained structure, after the
        relaxation if `relax_job` is supplied.
    strain_kwargs
        Additional keyword arguments to pass to
        [quacc.atoms.deformation.make_volume_strains_from_bulk][]
    include_undeformed
        Whether to include the energy of `undeformed_result` as a data point in the
        fit. Set this to False if `undeformed_result` is not directly comparable to
        the volume-strained calculations (e.g. the final step of a cell relaxation
        when the strained structures are run with a static calculation).

    Returns
    -------
    EOSSchema
        Results of the undeformed and volume-strained calculations and the
        equation of state fit. The volume-strained results are grouped under
        "relax" and "static", and a list is empty when its job is not supplied.
        See the return type-hint for the data structure.
    """
    if relax_job is None and static_job is None:
        raise ValueError("At least one of relax_job or static_job must be supplied.")

    strain_kwargs = strain_kwargs or {}
    strained_atoms = make_volume_strains_from_bulk(
        undeformed_result["atoms"], **strain_kwargs
    )
    if len(strained_atoms) + include_undeformed < 4:
        raise ValueError(
            "At least 4 energy-volume points are needed to fit an equation of state."
        )

    deformed_results: dict[str, list[Schema]] = {"relax": [], "static": []}
    for atoms in strained_atoms:
        final_atoms = atoms
        if relax_job is not None:
            result = relax_job(atoms, relax_cell=False)
            deformed_results["relax"].append(result)
            final_atoms = result["atoms"]

        if static_job is not None:
            result = static_job(final_atoms)
            deformed_results["static"].append(result)

    fit_results: list[RunSchema] = [
        *deformed_results["static" if static_job is not None else "relax"]
    ]
    if include_undeformed:
        fit_results.insert(0, undeformed_result)

    return {
        "undeformed_result": undeformed_result,
        "deformed_results": deformed_results,
        "eos_doc": _fit_eos(fit_results),
    }


@job
def _fit_eos(results: list[RunSchema]) -> EOSDoc:
    """
    Fit energy-volume data to several equations of state.

    Parameters
    ----------
    results
        Results of the calculations to fit, each with an energy and a structure.

    Returns
    -------
    EOSDoc
        An emmet (i.e. MaterialsProject) EOSDoc with a fit for each supported
        equation of state, including the equilibrium volume (Å³) and bulk
        modulus (GPa).
    """
    results = sorted(results, key=lambda result: result["atoms"].get_volume())
    return EOSDoc.from_ev_data(
        volumes=[result["atoms"].get_volume() for result in results],
        energies=[result["results"]["energy"] for result in results],
    )

"""
MOF-off-compatible VASP static calculation recipe.
"""

from __future__ import annotations

from importlib.util import find_spec
from typing import TYPE_CHECKING

from ase.units import Bohr
from monty.dev import requires

from quacc import flow, job
from quacc.recipes.vasp.matpes import matpes_static_job
from quacc.utils.dicts import recursive_dict_merge
from quacc.wflow_tools.customizers import customize_funcs

has_atomate2 = bool(find_spec("atomate2"))

if TYPE_CHECKING:
    from collections.abc import Callable
    from typing import Any, Literal, TypedDict

    from ase.atoms import Atoms

    from quacc.types import SourceDirectory, VaspSchema

    class MOFOffStaticFlowSchema(TypedDict):
        """Type hint associated with the MOF-off static flow."""

        pbe: VaspSchema
        r2scan: VaspSchema


@job
@requires(has_atomate2, "atomate2 is not installed. Run `pip install quacc[mp]`")
def mof_off_static_job(
    atoms: Atoms,
    *,
    level: Literal["PBE", "r2SCAN"],
    dispersion: Literal["D3BJ", "D4"] | None = None,
    prev_dir: SourceDirectory | None = None,
    **calc_kwargs: Any,
) -> VaspSchema:
    """
    Function to run a MOF-Off-compatible static calculation.

    Parameters
    ----------
    atoms
        Atoms object
    level
        The level of theory: "PBE", "r2SCAN"
    dispersion
        Dispersion corrections: None, "D3BJ", "D4"
    prev_dir
        A previous directory for a prior step in the workflow.
    **calc_kwargs
        Custom kwargs for the Vasp calculator. Set a value to
        `None` to remove a pre-existing key entirely. For a list of available
        keys, refer to [quacc.calculators.vasp.vasp.Vasp][]. All of the ASE
        Vasp calculator keyword arguments are supported.

    Returns
    -------
    VaspSchema
        Dictionary of results from [quacc.schemas.vasp.VaspSummarize.run][].
        See the type-hint for the data structure.
    """
    default_parameters = {
        "kspacing": 0.4,
        "incar_copilot": "critical",
        "use_improvements": True,
        "write_extra_files": True,
    }
    if dispersion:
        disp = dispersion.lower()
        if disp == "d3bj":
            default_parameters |= {
                "ivdw": 15,
                "sdftd3_damping": "rational",
                "vdw_radius": 95 * Bohr,
            }
        elif disp == "d4":
            default_parameters["ivdw"] = 13
    calc_flags = recursive_dict_merge(default_parameters, calc_kwargs)

    return matpes_static_job(
        atoms, level=level, auto_ispin=True, prev_dir=prev_dir, **calc_flags
    )


@flow
@requires(has_atomate2, "atomate2 is not installed. Run `pip install quacc[mp]`")
def mof_off_static_flow(
    atoms: Atoms,
    dispersion: Literal["D3BJ", "D4"] | None = None,
    job_params: dict[str, dict[str, Any]] | None = None,
    job_decorators: dict[str, Callable | None] | None = None,
) -> MOFOffStaticFlowSchema:
    """Run consecutive MOF-off-compatible PBE and r2SCAN static calculations.

    Parameters
    ----------
    atoms
        Atoms object.
    dispersion
        Dispersion correction to apply to both static calculations: None,
        "D3BJ", or "D4".
    job_params
        Custom parameters for the ``pbe_static_job`` and ``r2scan_static_job``
        steps. Use the ``all`` key to customize both steps.
    job_decorators
        Custom decorators for the ``pbe_static_job`` and ``r2scan_static_job``
        steps. Use the ``all`` key to customize both steps.

    Returns
    -------
    MOFOffStaticFlowSchema
        Dictionary containing the PBE and r2SCAN results.
    """
    pbe_static_job, r2scan_static_job = customize_funcs(
        ["pbe_static_job", "r2scan_static_job"],
        [mof_off_static_job, mof_off_static_job],
        param_defaults={
            "pbe_static_job": {"level": "PBE", "dispersion": dispersion},
            "r2scan_static_job": {"level": "r2SCAN", "dispersion": dispersion},
        },
        param_swaps=job_params,
        decorators=job_decorators,
    )

    pbe_results = pbe_static_job(atoms)
    r2scan_results = r2scan_static_job(
        pbe_results["atoms"], prev_dir=pbe_results["dir_name"]
    )

    return {"pbe": pbe_results, "r2scan": r2scan_results}

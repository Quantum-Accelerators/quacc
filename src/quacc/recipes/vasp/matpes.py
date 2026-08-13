"""
MatPES-compatible recipes

!!! Important

    Make sure that you use the MatPES-compatible pseudpotential
    versions (i.e. v.64)
"""

from __future__ import annotations

from importlib.util import find_spec
from typing import TYPE_CHECKING

from monty.dev import requires

from quacc import flow, job
from quacc.calculators.vasp.params import MPtoASEConverter
from quacc.recipes.vasp._base import run_and_summarize
from quacc.wflow_tools.customizers import customize_jobs

has_atomate2 = bool(find_spec("atomate2"))

if TYPE_CHECKING:
    from collections.abc import Callable
    from typing import Any, Literal, TypedDict

    from ase.atoms import Atoms

    from quacc.types import SourceDirectory, VaspSchema

    class MatPESStaticFlowSchema(TypedDict):
        """Type hint associated with the MatPES static flow."""

        pbe: VaspSchema
        r2scan: VaspSchema

@job
@requires(has_atomate2, "atomate2 is not installed. Run `pip install quacc[mp]`")
def matpes_static_job(
    atoms: Atoms,
    *,
    level: Literal["PBE", "r2SCAN", "HSE06"],
    kspacing: float | None = 0.22,
    use_improvements: bool = False,
    write_extra_files: bool = False,
    auto_ispin: bool = False,
    prev_dir: SourceDirectory | None = None,
    **calc_kwargs: Any,
) -> VaspSchema:
    """
    Function to run a MatPES-compatible static calculation.

    Parameters
    ----------
    atoms
        Atoms object
    level
        The level of theory: "PBE", "r2SCAN", "HSE06"
    kspacing
        The KSPACING parameter to use. Default: 0.22 as in the MatPES
        paper. This is likely too expensive in many cases.
    use_improvements
        Whether to make the following improvements to the VASP settings:
        ALGO = All, EFERMI = MIDGAP, GGA_COMPAT = False, ISEARCH = 1,
        and ENAUG deleted.
    write_extra_files
        Whether to write out the following IO files: LELF = True and NEDOS = 3001.
    auto_ispin
        If generating input set from a previous calculation, this controls whether
        to disable magnetisation (ISPIN = 1) if the absolute value of all magnetic
        moments are less than 0.02.
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
    from atomate2.vasp.jobs.matpes import MatPesGGAStaticMaker

    maker = MatPesGGAStaticMaker()
    maker.input_set_generator.auto_ispin = auto_ispin
    calc_defaults = MPtoASEConverter(atoms=atoms, prev_dir=prev_dir).convert_maker(
        maker
    )

    # Set the user-defined KSPACING
    calc_defaults |= {"kspacing": kspacing, "incar_copilot": "critical"}

    # Set some parameters that we think are improvements to MatPES
    if use_improvements:
        calc_defaults |= {
            "algo": "all",
            "efermi": "midgap",
            "enaug": None,
            "gga_compat": False,
            "isearch": 1,
        }

    # Write out optional files
    if write_extra_files:
        calc_defaults |= {"lelf": True, "nedos": 3001}

    # Set the level of theory
    del calc_defaults["gga"]
    if level.lower() == "pbe":
        calc_defaults |= {"xc": "pbe", "lwave": True}
    elif level.lower() == "r2scan":
        calc_defaults |= {"xc": "r2scan"}
    elif level.lower() == "hse06":
        calc_defaults |= {"algo": "normal", "xc": "hse06"}
        calc_defaults.pop("isearch", None)
    else:
        raise ValueError(f"Unsupported value for {level}")

    return run_and_summarize(
        atoms,
        calc_defaults=calc_defaults,
        calc_swaps=calc_kwargs,
        additional_fields={"name": f"MatPES {level} Static"},
        copy_files=[{"source": prev_dir, "filenames": ["WAVECAR*"]}]
        if prev_dir
        else None,
    )


@flow
@requires(has_atomate2, "atomate2 is not installed. Run `pip install quacc[mp]`")
def matpes_static_flow(
    atoms: Atoms,
    job_params: dict[str, dict[str, Any]] | None = None,
    job_decorators: dict[str, Callable | None] | None = None,
) -> MatPESStaticFlowSchema:
    """Run consecutive MatPES-compatible PBE and r2SCAN static calculations.

    Parameters
    ----------
    atoms
        Atoms object.
    job_params
        Custom parameters for the ``pbe_static_job`` and ``r2scan_static_job``
        steps. Use the ``all`` key to customize both steps.
    job_decorators
        Custom decorators for the ``pbe_static_job`` and ``r2scan_static_job``
        steps. Use the ``all`` key to customize both steps.

    Returns
    -------
    MatPESStaticFlowSchema
        Dictionary containing the PBE and r2SCAN results.
    """
    pbe_static_job, r2scan_static_job = customize_jobs(
        {"pbe_static_job": matpes_static_job, "r2scan_static_job": matpes_static_job},
        param_defaults={
            "pbe_static_job": {"level": "PBE"},
            "r2scan_static_job": {"level": "r2SCAN"},
        },
        param_swaps=job_params,
        decorators=job_decorators,
    ).values()

    pbe_results = pbe_static_job(atoms)
    r2scan_results = r2scan_static_job(
        pbe_results["atoms"], prev_dir=pbe_results["dir_name"]
    )

    return {"pbe": pbe_results, "r2scan": r2scan_results}

@job
@requires(has_atomate2, "atomate2 is not installed. Run `pip install quacc[mp]`")
def fastpes_static_job(
    atoms: Atoms,
    *,
    hf_exchange: float = 0.0,
    prev_dir: SourceDirectory | None = None,
    **calc_kwargs: Any,
) -> VaspSchema:
    """
    Function to run a FastPES-compatible static calculation.

    Parameters
    ----------
    atoms
        Atoms object
    hf_exchange
        The fraction of Hartree-Fock exchange, from 0.0 (PBE) to 1.0.
        HSE06 is 0.25.
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
    return matpes_static_job(atoms, level="PBE" if hf_exchange == 0.0 else "HSE06",kspacing=0.40,use_improvements=True,encut=435, isif=2,setups=_fastpes_paw,aexx=hf_exchange,additional_fields={"name": f"FastPES {hf_exchange} Static"}"

_fastpes_paw = {'Ac': '', 'Ag': '_GW', 'Al': '', 'Am': '', 'Ar': '_GW', 'As': '_GW', 'At': '', 'Au': '_GW', 'B': '_GW', 'Ba': '_sv_GW', 'Be': '', 'Bi': '_GW', 'Br': '_GW', 'C': '_s_GW', 'Ca': '_pv', 'Cd': '', 'Ce': '_3', 'Cf': '', 'Cl': '_GW', 'Cm': '', 'Co': '_GW', 'Cr': '', 'Cs': '_sv', 'Cu': '', 'Dy': '_3', 'Er': '_2', 'Eu': '_2', 'F': '_s', 'Fe': '_GW', 'Fr': '_sv', 'Ga': '', 'Gd': '_3', 'Ge': '', 'H': '_GW', 'He': '_GW', 'Hf': '', 'Hg': '', 'Ho': '_3', 'I': '', 'In': '', 'Ir': '', 'K': '_pv', 'Kr': '_GW', 'La': '_s', 'Li': '', 'Lu': '_3', 'Mg': '', 'Mn': '_GW', 'Mo': '', 'N': '_s_GW', 'Na': '', 'Nb': '_pv', 'Nd': '_3', 'Ne': '_s_GW', 'Ni': '', 'Np': '', 'O': '_s_GW', 'Os': '', 'P': '_GW', 'Pa': '_s', 'Pb': '', 'Pd': '', 'Pm': '_3', 'Po': '', 'Pr': '_3', 'Pt': '_GW', 'Pu': '', 'Ra': '_sv', 'Rb': '_pv', 'Re': '', 'Rh': '_GW', 'Rn': '', 'Ru': '', 'S': '_GW', 'Sb': '', 'Sc': '', 'Se': '', 'Si': '_GW', 'Sm': '_3', 'Sn': '', 'Sr': '_sv', 'Ta': '', 'Tb': '_3', 'Tc': '', 'Te': '', 'Th': '_s', 'Ti': '', 'Tl': '', 'Tm': '_3', 'U': '', 'V': '', 'W': '', 'Xe': '_GW', 'Y': '_sv', 'Yb': '_2', 'Zn': '_GW', 'Zr': '_sv'}

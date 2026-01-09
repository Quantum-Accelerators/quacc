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

from quacc import job
from quacc.calculators.vasp.params import MPtoASEConverter
from quacc.recipes.vasp._base import run_and_summarize
from quacc.wflow_tools.job_argument import Copy

has_atomate2 = bool(find_spec("atomate2"))

if TYPE_CHECKING:
    from typing import Literal

    from ase.atoms import Atoms

    from quacc.types import SourceDirectory, VaspSchema


@job
@requires(has_atomate2, "atomate2 is not installed. Run `pip install quacc[mp]`")
def matpes_static_job(
    atoms: Atoms,
    *,
    level: Literal["PBE", "r2SCAN", "HSE06"],
    kspacing: float | None = 0.22,
    use_improvements: bool = True,
    write_extra_files: bool = True,
    prev_dir: SourceDirectory | None = None,
    **calc_kwargs,
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
    maker.input_set_generator.auto_ispin = True
    calc_defaults = MPtoASEConverter(atoms=atoms, prev_dir=prev_dir).convert_maker(
        maker
    )

    # Set the user-defined KSPACING
    calc_defaults["kspacing"] = kspacing

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
        copy_files=Copy(src_files=prev_dir, files=["WAVECAR*"]) if prev_dir else None,
    )

"""
MOF-off-compatible VASP static calculation recipe.
"""

from __future__ import annotations

from importlib.util import find_spec
from typing import TYPE_CHECKING

from monty.dev import requires

from quacc import job
from quacc.recipes.vasp.matpes import matpes_static_job
from quacc.utils.dicts import recursive_dict_merge

has_atomate2 = bool(find_spec("atomate2"))

if TYPE_CHECKING:
    from typing import Literal

    from ase.atoms import Atoms

    from quacc.types import SourceDirectory, VaspSchema


@job
@requires(has_atomate2, "atomate2 is not installed. Run `pip install quacc[mp]`")
def mof_off_static_job(
    atoms: Atoms,
    *,
    level: Literal["PBE", "r2SCAN"],
    prev_dir: SourceDirectory | None = None,
    **calc_kwargs,
) -> VaspSchema:
    """
    Function to run a MOF-Off-compatible static calculation.

    Parameters
    ----------
    atoms
        Atoms object
    level
        The level of theory: "PBE", "r2SCAN"
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
        "use_improvements": True,
        "write_extra_files": True,
    }
    calc_flags = recursive_dict_merge(default_parameters, calc_kwargs)

    return matpes_static_job(atoms, level=level, prev_dir=prev_dir, **calc_flags)

"""
Custom types used throughout quacc.
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING, Any, Callable, Literal, TypedDict, Union

if TYPE_CHECKING:
    from numpy.typing import ArrayLike, NDArray

Filenames = Union[str, Path, list[Union[str, Path]]]
SourceDirectory = Union[str, Path]


class DefaultSetting:
    """
    Type hint for when a default setting will be applied
    """


class PmgKpts(TypedDict, total=False):
    """
    Type hint for `pmg_kpts` in [quacc.utils.kpts.convert_pmg_kpts][].
    """

    line_density: float
    kppvol: float
    kppa: float
    length_densities: tuple[float, float, float]


class AdsSiteFinderKwargs(TypedDict, total=False):
    """
    Type hint for `ads_site_finder_kwargs` in [quacc.atoms.slabs.make_adsorbate_structures][].
    """

    selective_dynamics: bool  # default = False
    height: float  # default = 0.9
    mi_vec: ArrayLike | None  # default = None


class FindAdsSitesKwargs(TypedDict, total=False):
    """
    Type hint for `find_ads_sites_kwargs` in [quacc.atoms.slabs.make_adsorbate_structures][].
    """

    distance: float  # default = 2.0
    put_inside: bool  # default = True
    symm_reduce: float  # default = 1e-2
    near_reduce: float  # default = 1e-2
    positions: list[
        Literal["ontop", "bridge", "hollow", "subsurface"]
    ]  # default: ["ontop", "bridge", "hollow"]
    no_obtuse_hollow: bool  # default = True


class QchemResults(TypedDict, total=False):
    """
    Type hint for the `results` attribute in [quacc.calculators.qchem.qchem.QChem][].
    """

    energy: float  # electronic energy in eV
    forces: NDArray  # forces in eV/A
    hessian: NDArray  # Hessian in eV/A^2/amu
    taskdoc: dict[str, Any]  # Output from `emmet.core.qc_tasks.TaskDoc`


class VaspJobKwargs(TypedDict, total=False):
    """
    Type hint for `vasp_job_kwargs` in in [quacc.calculators.vasp.vasp_custodian.run_custodian][].
    """

    output_file: str  # default = "vasp.out"
    stderr_file: str  # default = "std_err.txt"
    suffix: str  # default = ""
    final: bool  # default = True
    backup: bool  # default = True
    auto_npar: bool  # default = False
    auto_gamma: bool  # default = True
    settings_override: dict | None  # default = None
    copy_magmom: bool  # default = False
    auto_continue: bool  # default = False


class VaspCustodianKwargs(TypedDict, total=False):
    """
    Type hint for `custodian_kwargs` in [quacc.calculators.vasp.vasp_custodian.run_custodian][].
    """

    max_errors_per_job: int | None  # default = None
    polling_time_step: int  # default = 10
    monitor_freq: int  # default = 10
    skip_over_errors: bool  # default = False
    gzipped_output: bool  # default = False
    checkpoint: bool  # default = False
    terminate_func: Callable | None  # default = None
    terminate_on_nonzero_returncode: bool  # default = False

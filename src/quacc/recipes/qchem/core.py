"""Core recipes for Q-Chem."""

from __future__ import annotations

from importlib.util import find_spec
from typing import TYPE_CHECKING

from quacc import job
from quacc.recipes.qchem._base import run_and_summarize, run_and_summarize_opt
from quacc.utils.dicts import recursive_dict_merge

has_sella = bool(find_spec("sella"))

if has_sella:
    from sella import Sella

if TYPE_CHECKING:
    from typing import Any

    from ase.atoms import Atoms

    from quacc.types import Filenames, OptParams, OptSchema, RunSchema, SourceDirectory

_BASE_SET = {
    "rem": {
        "gen_scfman": True,
        "xc_grid": 3,
        "thresh": 14,
        "s2thresh": 16,
        "scf_algorithm": "diis",
        "resp_charges": True,
        "symmetry": False,
        "sym_ignore": True,
    }
}


@job
def static_job(
    atoms: Atoms,
    charge: int = 0,
    spin_multiplicity: int = 1,
    method: str | None = "wb97mv",
    basis: str | None = "def2-tzvpd",
    copy_files: SourceDirectory | dict[SourceDirectory, Filenames] | None = None,
    additional_fields: dict[str, Any] | None = None,
    **calc_kwargs,
) -> RunSchema:
    """
    Carry out a single-point calculation.

    Parameters
    ----------
    atoms
        Atoms object
    charge
        Total charge of the system.
    spin_multiplicity
        Multiplicity of the system.
    method
        DFT exchange-correlation functional or other electronic structure
        method.
    basis
        Basis set.
    copy_files
        Files to copy (and decompress) from source to the runtime directory.
    additional_fields
        Additional fields to add to the results dictionary.
    **calc_kwargs
        Custom kwargs for the calculator. Set a value to `quacc.Remove` to remove
        a pre-existing key entirely. See [quacc.calculators.qchem.qchem.QChem][] for more
        details.

    Returns
    -------
    RunSchema
        Dictionary of results from [quacc.schemas.ase.Summarize.run][].
        See the type-hint for the data structure.
    """
    calc_defaults = recursive_dict_merge(
        _BASE_SET, {"rem": {"job_type": "force", "method": method, "basis": basis}}
    )

    return run_and_summarize(
        atoms,
        charge=charge,
        spin_multiplicity=spin_multiplicity,
        calc_defaults=calc_defaults,
        calc_swaps=calc_kwargs,
        additional_fields={"name": "Q-Chem Static"} | (additional_fields or {}),
        copy_files=copy_files,
    )


@job
def relax_job(
    atoms: Atoms,
    charge: int = 0,
    spin_multiplicity: int = 1,
    method: str = "wb97mv",
    basis: str = "def2-svpd",
    opt_params: OptParams | None = None,
    copy_files: SourceDirectory | dict[SourceDirectory, Filenames] | None = None,
    additional_fields: dict[str, Any] | None = None,
    **calc_kwargs,
) -> OptSchema:
    """
    Optimize aka "relax" a molecular structure with an ASE optimizer.

    Parameters
    ----------
    atoms
        Atoms object
    charge
        Charge of the system.
    spin_multiplicity
        Multiplicity of the system.
    method
        DFT exchange-correlation functional or other electronic structure
        method.
    basis
        Basis set.
    opt_params
        Dictionary of custom kwargs for the optimization process. For a list
        of available keys, refer to [quacc.runners.ase.Runner.run_opt][].
    copy_files
        Files to copy (and decompress) from source to the runtime directory.
    additional_fields
        Additional fields to add to the results dictionary.
    **calc_kwargs
        Custom kwargs for the calculator. Set a value to `quacc.Remove` to remove
        a pre-existing key entirely. See [quacc.calculators.qchem.qchem.QChem][] for more
        details.

    Returns
    -------
    OptSchema
        Dictionary of results from [quacc.schemas.ase.Summarize.opt][].
        See the type-hint for the data structure.
    """
    calc_defaults = recursive_dict_merge(
        _BASE_SET, {"rem": {"job_type": "force", "method": method, "basis": basis}}
    )
    opt_defaults = {"optimizer": Sella} if has_sella else {}

    return run_and_summarize_opt(
        atoms,
        charge=charge,
        spin_multiplicity=spin_multiplicity,
        calc_defaults=calc_defaults,
        calc_swaps=calc_kwargs,
        opt_defaults=opt_defaults,
        opt_params=opt_params,
        additional_fields={"name": "Q-Chem Relax"} | (additional_fields or {}),
        copy_files=copy_files,
    )


@job
def freq_job(
    atoms: Atoms,
    charge: int = 0,
    spin_multiplicity: int = 1,
    method: str = "wb97mv",
    basis: str = "def2-svpd",
    copy_files: SourceDirectory | dict[SourceDirectory, Filenames] | None = None,
    additional_fields: dict[str, Any] | None = None,
    **calc_kwargs,
) -> RunSchema:
    """
    Perform a frequency calculation on a molecular structure.

    Parameters
    ----------
    atoms
        Atoms object
    charge
        Charge of the system.
    spin_multiplicity
        Multiplicity of the system.
    method
        DFT exchange-correlation functional or other electronic structure
        method.
    basis
        Basis set.
    copy_files
        Files to copy (and decompress) from source to the runtime directory.    **calc_kwargs
        Custom kwargs for the calculator. Set a value to `quacc.Remove` to remove
        a pre-existing key entirely. See [quacc.calculators.qchem.qchem.QChem][] for more
        details.
    additional_fields
        Additional fields to add to the results dictionary.

    Returns
    -------
    RunSchema
        Dictionary of results from [quacc.schemas.ase.Summarize.run][]
    """
    calc_defaults = recursive_dict_merge(
        _BASE_SET, {"rem": {"job_type": "freq", "method": method, "basis": basis}}
    )

    return run_and_summarize(
        atoms,
        charge=charge,
        spin_multiplicity=spin_multiplicity,
        calc_defaults=calc_defaults,
        calc_swaps=calc_kwargs,
        copy_files=copy_files,
        additional_fields={"name": "Q-Chem Frequency"} | (additional_fields or {}),
    )

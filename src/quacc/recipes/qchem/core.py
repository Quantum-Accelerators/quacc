"""Core recipes for Q-Chem."""
from __future__ import annotations

from typing import TYPE_CHECKING

from ase.optimize import FIRE

from quacc import job
from quacc.recipes.qchem._base import base_fn, base_opt_fn
from quacc.utils.dicts import merge_dicts

try:
    from sella import Sella

    has_sella = True
except ImportError:
    has_sella = False

if TYPE_CHECKING:
    from pathlib import Path
    from typing import Any

    from ase.atoms import Atoms

    from quacc.schemas._aliases.ase import OptSchema, RunSchema

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
    copy_files: str | Path | list[str | Path] | None = None,
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
        File(s) to copy to the runtime directory. If a directory is provided, it will be recursively unpacked.
    **calc_kwargs
        Custom kwargs for the calculator. Set a value to `None` to remove
        a pre-existing key entirely. See [quacc.calculators._qchem_legacy.qchem.QChem][] for more
        details.

    Returns
    -------
    RunSchema
        Dictionary of results from [quacc.schemas.ase.summarize_run][]
    """

    calc_defaults = merge_dicts(
        _BASE_SET, {"rem": {"job_type": "force", "method": method, "basis": basis}}
    )

    return base_fn(
        atoms,
        charge=charge,
        spin_multiplicity=spin_multiplicity,
        calc_defaults=calc_defaults,
        calc_swaps=calc_kwargs,
        additional_fields={"name": "Q-Chem Static"},
        copy_files=copy_files,
    )


@job
def relax_job(
    atoms: Atoms,
    charge: int = 0,
    spin_multiplicity: int = 1,
    method: str = "wb97mv",
    basis: str = "def2-svpd",
    opt_params: dict[str, Any] | None = None,
    copy_files: str | Path | list[str | Path] | None = None,
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
        method. Defaults to wB97M-V.
    basis
        Basis set. Defaults to def2-SVPD.
    opt_params
        Dictionary of custom kwargs for the optimization process. Set a value
        to `None` to remove a pre-existing key entirely. For a list of available
        keys, refer to [quacc.runners.ase.run_opt][].
    copy_files
        File(s) to copy to the runtime directory. If a directory is provided, it will be recursively unpacked.
    **calc_kwargs
        Custom kwargs for the calculator. Set a value to `None` to remove
        a pre-existing key entirely. See [quacc.calculators._qchem_legacy.qchem.QChem][] for more
        details.

    Returns
    -------
    OptSchema
        Dictionary of results from [quacc.schemas.ase.summarize_opt_run][]
    """

    calc_defaults = merge_dicts(
        _BASE_SET, {"rem": {"job_type": "force", "method": method, "basis": basis}}
    )
    opt_defaults = {
        "fmax": 0.01,
        "max_steps": 1000,
        "optimizer": Sella if has_sella else FIRE,
    }

    return base_opt_fn(
        atoms,
        charge=charge,
        spin_multiplicity=spin_multiplicity,
        calc_defaults=calc_defaults,
        calc_swaps=calc_kwargs,
        opt_defaults=opt_defaults,
        opt_params=opt_params,
        additional_fields={"name": "Q-Chem Optimization"},
        copy_files=copy_files,
    )


@job
def freq_job(
    atoms: Atoms,
    charge: int = 0,
    spin_multiplicity: int = 1,
    method: str = "wb97mv",
    basis: str = "def2-svpd",
    copy_files: str | Path | list[str | Path] | None = None,
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
        method. Defaults to wB97M-V.
    basis
        Basis set. Defaults to def2-SVPD.
    copy_files
        File(s) to copy to the runtime directory. If a directory is provided, it will be recursively unpacked.
    **calc_kwargs
        Custom kwargs for the calculator. Set a value to `None` to remove
        a pre-existing key entirely. See [quacc.calculators._qchem_legacy.qchem.QChem][] for more
        details.

    Returns
    -------
    RunSchema
        Dictionary of results from [quacc.schemas.ase.summarize_run][]
    """

    calc_defaults = merge_dicts(
        _BASE_SET, {"rem": {"job_type": "freq", "method": method, "basis": basis}}
    )

    return base_fn(
        atoms,
        charge=charge,
        spin_multiplicity=spin_multiplicity,
        calc_defaults=calc_defaults,
        calc_swaps=calc_kwargs,
        copy_files=copy_files,
        additional_fields={"name": "Q-Chem Frequency"},
    )

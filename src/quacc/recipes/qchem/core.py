"""Core recipes for the Q-Chem."""
from __future__ import annotations

from typing import TYPE_CHECKING

from ase.optimize import FIRE

from quacc import job
from quacc.calculators.qchem import QChem
from quacc.runners.ase import run_calc, run_opt
from quacc.schemas.ase import summarize_opt_run, summarize_run
from quacc.utils.dicts import merge_dicts

try:
    from sella import Sella

    has_sella = True
except ImportError:
    has_sella = False

if TYPE_CHECKING:
    from typing import Any

    from ase import Atoms

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
    copy_files: list[str] | None = None,
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
        Files to copy to the runtime directory.
    **calc_kwargs
        Custom kwargs for the calculator. Set a value to `None` to remove
        a pre-existing key entirely. See [quacc.calculators.qchem.qchem.QChem][] for more
        details.

        !!! Info "Calculator defaults"

            ```python
            {
                "rem": {
                    "job_type": "force",
                    "method": method,
                    "basis": basis,
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
            ```

    Returns
    -------
    RunSchema
        Dictionary of results from [quacc.schemas.ase.summarize_run][]
    """
    defaults = _BASE_SET | {
        "job_type": "force",
        "method": method,
        "basis": basis,
    }

    return _base_job(
        atoms,
        charge,
        spin_multiplicity,
        defaults=defaults,
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
    copy_files: list[str] | None = None,
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

        !!! Info "Optimizer defaults"

            ```python
            {"fmax": 0.01, "max_steps": 1000, "optimizer": Sella if has_sella else FIRE}
            ```
    copy_files
        Files to copy to the runtime directory.
    **calc_kwargs
        Custom kwargs for the calculator. Set a value to `None` to remove
        a pre-existing key entirely. See [quacc.calculators.qchem.qchem.QChem][] for more
        details.

        !!! Info "Calculator defaults"

            ```python
            _BASE_SET | {
                    "job_type": "force",
                    "method": method,
                    "basis": basis,
                }
            ```

    Returns
    -------
    OptSchema
        Dictionary of results from [quacc.schemas.ase.summarize_opt_run][]
    """

    calc_defaults = _BASE_SET | {
        "job_type": "force",
        "method": method,
        "basis": basis,
    }

    opt_defaults = {
        "fmax": 0.01,
        "max_steps": 1000,
        "optimizer": Sella if has_sella else FIRE,
    }

    return _base_opt_job(
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
    copy_files: list[str] | None = None,
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
        Files to copy to the runtime directory.
    **calc_kwargs
        Custom kwargs for the calculator. Set a value to `None` to remove
        a pre-existing key entirely. See [quacc.calculators.qchem.qchem.QChem][] for more
        details.

        !!! Info "Calculator defaults"

            ```python
            _BASE_SET | {
                    "job_type": "freq",
                    "method": method,
                    "basis": basis,
                }
            ```

    Returns
    -------
    RunSchema
        Dictionary of results from [quacc.schemas.ase.summarize_run][]
    """

    defaults = _BASE_SET | {
        "job_type": "freq",
        "method": method,
        "basis": basis,
    }

    return _base_job(
        atoms,
        charge=charge,
        spin_multiplicity=spin_multiplicity,
        defaults=defaults,
        calc_swaps=calc_kwargs,
        copy_files=copy_files,
        additional_fields={"name": "Q-Chem Frequency"},
    )


def _base_job(
    atoms: Atoms,
    charge: int = 0,
    spin_multiplicity: int = 1,
    defaults: dict[str, Any] | None = None,
    calc_swaps: dict[str, Any] | None = None,
    additional_fields: dict[str, Any] | None = None,
    copy_files: list[str] | None = None,
) -> RunSchema:
    """
    Base job function used for Q-Chem recipes that don't rely on ASE optimizers or other
    ASE dynamics classes.

    Parameters
    ----------
    atoms
        Atoms object
    charge
        Charge of the system.
    spin_multiplicity
        Multiplicity of the system.
    defaults
        The default parameters for the recipe.
    calc_swaps
        Dictionary of custom kwargs for the calculator. Set a value to `None` to remove
        a pre-existing key entirely. See [quacc.calculators.qchem.qchem.QChem][] for more
        details.
    additional_fields
        Any additional fields to set in the summary.
    copy_files
        Files to copy to the runtime directory.

    Returns
    -------
    RunSchema
        Dictionary of results from [quacc.schemas.ase.summarize_run][]
    """

    qchem_flags = merge_dicts(defaults, calc_swaps)

    atoms.calc = QChem(atoms, **qchem_flags)
    final_atoms = run_calc(atoms, copy_files=copy_files)

    return summarize_run(
        final_atoms,
        input_atoms=atoms,
        charge_and_multiplicity=(charge, spin_multiplicity),
        additional_fields=additional_fields,
    )


def _base_opt_job(
    atoms: Atoms,
    charge: int = 0,
    spin_multiplicity: int = 1,
    calc_defaults: dict[str, Any] | None = None,
    calc_swaps: dict[str, Any] | None = None,
    opt_defaults: dict[str, Any] | None = None,
    opt_params: dict[str, Any] | None = None,
    additional_fields: dict[str, Any] | None = None,
    copy_files: list[str] | None = None,
) -> OptSchema:
    """
    Base function for Q-Chem recipes that involve ASE optimizers.

    Parameters
    ----------
    atoms
        Atoms object
    charge
        Charge of the system.
    spin_multiplicity
        Multiplicity of the system.
    calc_defaults
        Default arguments for the Q-Chem calculator.
    calc_swaps
        Dictionary of custom kwargs for the calculator. Set a value to `None` to remove
        a pre-existing key entirely. See [quacc.calculators.qchem.qchem.QChem][] for more
        details.
    opt_defaults
        Default arguments for the ASE optimizer.
    opt_params
        Dictionary of custom kwargs for [quacc.runners.ase.run_opt][]
    copy_files
        Files to copy to the runtime directory.

    Returns
    -------
    OptSchema
        Dictionary of results from [quacc.schemas.ase.summarize_opt_run][]
    """
    # TODO:
    #   - passing initial Hessian?

    qchem_flags = merge_dicts(calc_defaults, calc_swaps)
    opt_flags = merge_dicts(opt_defaults, opt_params)

    atoms.calc = QChem(atoms, **qchem_flags)
    dyn = run_opt(atoms, copy_files=copy_files, **opt_flags)

    return summarize_opt_run(
        dyn,
        charge_and_multiplicity=(charge, spin_multiplicity),
        additional_fields=additional_fields,
    )

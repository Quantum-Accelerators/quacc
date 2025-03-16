"""Core recipes for Gaussian."""

from __future__ import annotations

from typing import TYPE_CHECKING

from ase.calculators.gaussian import Gaussian

from quacc import job
from quacc.recipes.common.core import Recipe
from quacc.recipes.gaussian._util import create_gaussian_defaults

if TYPE_CHECKING:
    from typing import Any

    from ase.atoms import Atoms

    from quacc.types import Filenames, RunSchema, SourceDirectory


@job
def static_job(
    atoms: Atoms,
    charge: int = 0,
    spin_multiplicity: int = 1,
    xc: str = "wb97xd",
    basis: str = "def2tzvp",
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
        Charge of the system.
    spin_multiplicity
        Multiplicity of the system.
    xc
        Exchange-correlation functional
    basis
        Basis set
    copy_files
        Files to copy (and decompress) from source to the runtime directory.
    additional_fields
        Additional fields to add to the results dictionary.
    **calc_kwargs
        Custom kwargs for the Gaussian calculator. Set a value to
        `quacc.Remove` to remove a pre-existing key entirely. For a list of available
        keys, refer to the [ase.calculators.gaussian.Gaussian][] calculator.

    Returns
    -------
    RunSchema
        Dictionary of results
    """
    calc_defaults = create_gaussian_defaults(
        xc=xc, basis=basis, charge=charge, spin_multiplicity=spin_multiplicity
    )
    calc_defaults.update({"force": "", "gfinput": "", "ioplist": ["6/7=3", "2/9=2000"]})
    recipe = Recipe(Gaussian, calc_defaults)
    return recipe.calculate(
        atoms, copy_files=copy_files, additional_fields=additional_fields, **calc_kwargs
    )


@job
def relax_job(
    atoms: Atoms,
    charge: int,
    spin_multiplicity: int,
    xc: str = "wb97xd",
    basis: str = "def2tzvp",
    freq: bool = False,
    copy_files: SourceDirectory | dict[SourceDirectory, Filenames] | None = None,
    additional_fields: dict[str, Any] | None = None,
    **calc_kwargs,
) -> RunSchema:
    """
    Carry out a geometry optimization.

    Parameters
    ----------
    atoms
        Atoms object
    charge
        Charge of the system.
    spin_multiplicity
        Multiplicity of the system.
    xc
        Exchange-correlation functional
    basis
        Basis set
    freq
        If a frequency calculation should be carried out.
    copy_files
        Files to copy (and decompress) from source to the runtime directory.
    additional_fields
        Additional fields to add to the results
    **calc_kwargs
        Custom kwargs for the Gaussian calculator. Set a value to
        `quacc.Remove` to remove a pre-existing key entirely. For a list of available
        keys, refer to the [ase.calculators.gaussian.Gaussian][] calculator.

    Returns
    -------
    RunSchema
        Dictionary of results
    """
    calc_defaults = create_gaussian_defaults(
        xc=xc, basis=basis, charge=charge, spin_multiplicity=spin_multiplicity
    )
    calc_defaults |= {"opt": "", "ioplist": ["2/9=2000"]}  # ASE issue #660
    if freq:
        calc_defaults |= {"freq": ""}
    recipe = Recipe(Gaussian, calc_defaults)

    return recipe.calculate(
        atoms, copy_files=copy_files, additional_fields=additional_fields, **calc_kwargs
    )

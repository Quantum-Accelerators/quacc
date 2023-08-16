"""Core recipes for GULP"""
from __future__ import annotations

import warnings
from typing import TYPE_CHECKING

import covalent as ct
from ase.calculators.gulp import GULP

from quacc.schemas.ase import summarize_run
from quacc.schemas.atoms import fetch_atoms
from quacc.util.calc import run_calc
from quacc.util.dicts import get_parameters

if TYPE_CHECKING:
    from ase import Atoms

    from quacc.schemas.ase import RunSchema

GEOM_FILE_PBC = "gulp.cif"
GEOM_FILE_NOPBC = "gulp.xyz"


@ct.electron
def static_job(
    atoms: Atoms | dict,
    use_gfnff: bool = True,
    library: str | None = None,
    keyword_swaps: dict | None = None,
    option_swaps: dict | None = None,
    copy_files: list[str] | None = None,
) -> RunSchema:
    """
    Carry out a single-point calculation.
    Note: 'Conditions' are not yet natively supported.

    Parameters
    ----------
    atoms
        Atoms object or a dictionary with the key "atoms" and an Atoms object as the value
    use_gfnff
        True if (p)GFN-FF should be used; False if not.
    library
        Filename of the potential library file, if required.
    keyword_swaps
        Dictionary of custom keyword swap kwargs for the calculator.
    option_swaps
        Dictionary of custom option swap kwargs for the calculator.

    Returns
    -------
    RunSchema
        Dictionary of results from `quacc.schemas.ase.summarize_run`
    """
    atoms = fetch_atoms(atoms)

    default_keywords = {
        "gfnff": True if use_gfnff else None,
        "gwolf": True if use_gfnff and atoms.pbc.any() else None,
    }
    default_options = {
        "dump every gulp.res": True,
        f"output cif {GEOM_FILE_PBC}": True if atoms.pbc.any() else None,
        f"output xyz {GEOM_FILE_NOPBC}": None if atoms.pbc.any() else True,
    }

    keywords = get_parameters(default_keywords, swaps=keyword_swaps)
    options = get_parameters(default_options, swaps=option_swaps)

    gulp_keywords = " ".join(list(keywords.keys()))
    gulp_options = list(options.keys())

    atoms.calc = GULP(keywords=gulp_keywords, options=gulp_options, library=library)
    final_atoms = run_calc(
        atoms,
        geom_file=GEOM_FILE_PBC if atoms.pbc.any() else GEOM_FILE_NOPBC,
        copy_files=copy_files,
    )

    return summarize_run(
        final_atoms,
        input_atoms=atoms,
        additional_fields={"name": "GULP Static"},
    )


@ct.electron
def relax_job(
    atoms: Atoms | dict,
    use_gfnff: bool = True,
    library: str | None = None,
    relax_cell: bool = True,
    keyword_swaps: dict | None = None,
    option_swaps: dict | None = None,
    copy_files: list[str] | None = None,
) -> RunSchema:
    """
    Carry out a single-point calculation.
    Note: 'Conditions' are not yet natively supported.

    Parameters
    ----------
    atoms
        Atoms object or a dictionary with the key "atoms" and an Atoms object as the value
    use_gfnff
        True if (p)GFN-FF should be used; False if not.
    library
        Filename of the potential library file, if required.
    relax_cell
        True if the volume should be relaxed; False if not.
    keyword_swaps
        Dictionary of custom keyword swap kwargs for the calculator.
    option_swaps
        Dictionary of custom option swap kwargs for the calculator.
    copy_files
        Absolute paths to files to copy to the runtime directory.

    Returns
    -------
    dict
        Dictionary of results from `quacc.schemas.ase.summarize_run`
    """
    atoms = fetch_atoms(atoms)

    if relax_cell and not atoms.pbc.any():
        warnings.warn(
            "Volume relaxation requested but no PBCs found. Ignoring.", UserWarning
        )
        relax_cell = False

    default_keywords = {
        "opti": True,
        "gfnff": True if use_gfnff else None,
        "gwolf": True if use_gfnff and atoms.pbc.any() else None,
        "conp": True if relax_cell and atoms.pbc.any() else None,
        "conv": None if relax_cell and atoms.pbc.any() else True,
    }
    default_options = {
        "dump every gulp.res": True,
        f"output cif {GEOM_FILE_PBC}": True if atoms.pbc.any() else None,
        f"output xyz {GEOM_FILE_NOPBC}": None if atoms.pbc.any() else True,
    }

    keywords = get_parameters(default_keywords, swaps=keyword_swaps)
    options = get_parameters(default_options, swaps=option_swaps)

    gulp_keywords = " ".join(list(keywords.keys()))
    gulp_options = list(options.keys())

    atoms.calc = GULP(keywords=gulp_keywords, options=gulp_options, library=library)
    final_atoms = run_calc(
        atoms,
        geom_file=GEOM_FILE_PBC if atoms.pbc.any() else GEOM_FILE_NOPBC,
        copy_files=copy_files,
    )

    if not final_atoms.calc.get_opt_state():
        raise ValueError("Optimization did not converge!")

    return summarize_run(
        final_atoms,
        input_atoms=atoms,
        additional_fields={"name": "GULP Relax"},
    )

"""Utility functions for dealing with Atoms."""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING

import numpy as np

from quacc.atoms.core import copy_atoms, get_atoms_id

if TYPE_CHECKING:
    from ase.atoms import Atoms

logger = logging.getLogger(__name__)


def prep_next_run(atoms: Atoms, move_magmoms: bool = False) -> Atoms:
    """
    Prepares the Atoms object for a new run by stripping off the calculator and
    assigning a unique ID.

    Parameters
    ----------
    atoms
        Atoms object
    move_magmoms
        Whether to move the magnetic moments from the calculator results to the
        initial magnetic moments.

    Returns
    -------
    Atoms
        Updated Atoms object.
    """
    atoms = copy_atoms(atoms)

    if (
        move_magmoms
        and hasattr(atoms, "calc")
        and getattr(atoms.calc, "results", None) is not None
    ):
        atoms.set_initial_magnetic_moments(
            atoms.calc.results.get("magmoms", [0.0] * len(atoms))
        )

    # Clear off the calculator so we can run a new job. If we don't do this,
    # then something like atoms *= (2,2,2) still has a calculator attached,
    # which is a bit confusing.
    atoms.calc = None

    # Give the Atoms object a unique ID. This will be helpful for querying
    # later. Also store any old IDs somewhere else for future reference. Note:
    # Keep this at the end of the function so that the ID is assigned based on
    # the returned Atoms object.
    if atoms.info.get("_id", None) is not None:
        if atoms.info.get("_old_ids") is None:
            atoms.info["_old_ids"] = []
        atoms.info["_old_ids"].append(atoms.info["_id"])
    atoms.info["_id"] = get_atoms_id(atoms)

    return atoms


def set_magmoms(
    atoms: Atoms,
    elemental_mags_dict: dict[str, float] | None = None,
    elemental_mags_default: float = 1.0,
    copy_magmoms: bool = True,
    mag_cutoff: float | None = 0.05,
) -> Atoms:  # sourcery skip
    """
    Sets the initial magnetic moments in the Atoms object.

    This function deserves particular attention. The following logic is applied:
    - If there is a converged set of magnetic moments, those are moved to the
    initial magmoms if copy_magmoms is True.
    - If there is no converged set of magnetic moments but the user has set initial magmoms, those are simply used
    as is.
    - If there are no converged magnetic moments or initial magnetic
    moments, then the default magnetic moments from the preset
    elemental_mags_dict (if specified) are set as the initial magnetic moments.
    - For any of the above scenarios, if mag_cutoff is not None, the newly set
    initial magnetic moments are checked. If all have a magnitude below
    mag_cutoff, then they are all set to 0 (no spin polarization).

    Parameters
    ----------
    atoms
        Atoms object
    elemental_mags_dict
        Dictionary of elemental symbols and their corresponding magnetic moments
        to set. If None, no default values will be used.
    elemental_mags_default
        Default magnetic moment on an element if no magnetic moment is specified
        in the elemental_mags_dict. Only used if elemental_mags_dict is not
        None. This kwarg is mainly a convenience so that you don't need to list
        every single element in the elemental_mags_dict.
    copy_magmoms
        Whether to copy the magnetic moments from the converged set of magnetic
        moments to the initial magnetic moments.
    mag_cutoff
        Magnitude below which the magnetic moments are considered to be zero. If
        None, no cutoff will be applied

    Returns
    -------
    Atoms
        Atoms object
    """

    # Handle the magnetic moments Check if a prior job was run and pull the
    # prior magmoms
    if hasattr(atoms, "calc") and getattr(atoms.calc, "results", None) is not None:
        mags = atoms.calc.results.get("magmoms", [0.0] * len(atoms))
        # Note: It is important that we set mags to 0.0 here rather than None if
        # the calculator has no magmoms because: 1) ispin=1 might be set, and 2)
        # we do not want the preset magmoms to be used.
    else:
        mags = None

    # Check if the user has set any initial magmoms
    has_initial_mags = atoms.has("initial_magmoms")

    # If there are no initial magmoms set and this is not a follow-up job, we
    # may need to add some from the preset yaml.
    if mags is None:
        if not has_initial_mags:
            # If the preset dictionary has default magmoms, set those by
            # element. If the element isn't in the magmoms dict then set it to
            # mag_default.
            if elemental_mags_dict:
                initial_mags = np.array(
                    [
                        elemental_mags_dict.get(atom.symbol, elemental_mags_default)
                        for atom in atoms
                    ]
                )
                atoms.set_initial_magnetic_moments(initial_mags)
        else:
            pass
    elif copy_magmoms:
        atoms.set_initial_magnetic_moments(mags)

    # If all the set mags are below mag_cutoff, set them to 0
    if mag_cutoff:
        has_new_initial_mags = atoms.has("initial_magmoms")
        new_initial_mags = atoms.get_initial_magnetic_moments()
        if has_new_initial_mags and np.all(np.abs(new_initial_mags) < mag_cutoff):
            atoms.set_initial_magnetic_moments([0.0] * len(atoms))

    return atoms

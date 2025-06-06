"""Summarizer for phonopy."""

from __future__ import annotations

from importlib.util import find_spec
from typing import TYPE_CHECKING

from monty.dev import requires

from quacc import __version__
from quacc.schemas.atoms import atoms_to_metadata
from quacc.utils.dicts import clean_dict
from quacc.utils.files import get_uri

has_phonopy = bool(find_spec("phonopy"))

if TYPE_CHECKING:
    from pathlib import Path
    from typing import Any

    from ase.atoms import Atoms

    from quacc.types import PhononSchema

    if has_phonopy:
        from phonopy import Phonopy


@requires(has_phonopy, "This schema relies on phonopy")
def summarize_phonopy(
    phonon: Phonopy,
    input_atoms: Atoms,
    directory: str | Path,
    parameters: dict[str, Any] | None = None,
    additional_fields: dict[str, Any] | None = None,
) -> PhononSchema:
    """
    Summarize a Phonopy object.

    Parameters
    ----------
    phonon
        Phonopy object
    input_atoms
        Input atoms object
    directory
        Directory where the results are stored.
    parameters
        Calculator parameters used to generate the phonon object.
    additional_fields
        Additional fields to add to the document.

    Returns
    -------
    PhononSchema
        The PhononSchema.
    """
    additional_fields = additional_fields or {}

    inputs = {
        "parameters": parameters,
        "nid": get_uri(directory).split(":")[0],
        "dir_name": directory,
        "phonopy_metadata": {"version": phonon.version},
        "quacc_version": __version__,
    }

    results = {
        "results": {
            "thermal_properties": phonon.get_thermal_properties_dict(),
            "mesh_properties": phonon.get_mesh_dict(),
            "total_dos": phonon.get_total_dos_dict(),
            "force_constants": phonon.force_constants,
        }
    }

    atoms_metadata = atoms_to_metadata(input_atoms)
    unsorted_task_doc = atoms_metadata | inputs | results | additional_fields
    return clean_dict(unsorted_task_doc)

"""Summarizer for phonopy"""
from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

from quacc import SETTINGS
from quacc.schemas.atoms import atoms_to_metadata
from quacc.utils.dicts import clean_task_doc, recursive_dict_merge
from quacc.utils.files import get_uri
from quacc.wflow_tools.db import results_to_db

if TYPE_CHECKING:
    from typing import Any

    from ase.atoms import Atoms
    from maggma.core import Store
    from phonopy import Phonopy

    from quacc.schemas._aliases.phonons import PhononSchema


def summarize_phonopy(
    phonon: Phonopy,
    input_atoms: Atoms | None = None,
    parameters: dict[str, Any] | None = None,
    additional_fields: dict[str, Any] | None = None,
    store: Store | bool | None = None,
) -> PhononSchema:
    """
    Summarize a Phonopy object.

    Parameters
    ----------
    phonon
        Phonopy object
    input_atoms
        Input atoms object
    parameters
        Calculator parameters used to generate the phonon object.
    additional_fields
        Additional fields to add to the document.
    store
        Whether to store the document in the database.

    Returns
    -------
    PhononSchema
        The PhononSchema.
    """
    additional_fields = additional_fields or {}
    store = SETTINGS.PRIMARY_STORE if store is None else store

    uri = get_uri(Path.cwd())

    inputs = {
        "parameters": parameters,
        "phonopy_metadata": {"version": phonon.version},
        "nid": uri.split(":")[0],
        "dir_name": ":".join(uri.split(":")[1:]),
    }

    results = {
        "results": {
            "thermal_properties": phonon.get_thermal_properties_dict(),
            "mesh_properties": phonon.get_mesh_dict(),
            "force_constants": phonon.force_constants,
        }
    }
    phonon.save(settings={"force_constants": True})

    atoms_metadata = atoms_to_metadata(input_atoms) if input_atoms else {}
    unsorted_task_doc = recursive_dict_merge(
        atoms_metadata, inputs, results, additional_fields
    )
    task_doc = clean_task_doc(unsorted_task_doc)

    if store:
        results_to_db(store, task_doc)

    return task_doc

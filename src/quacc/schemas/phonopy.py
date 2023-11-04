"""Summarizer for phonopy"""
from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

from quacc import SETTINGS
from quacc.schemas.atoms import atoms_to_metadata
from quacc.utils.dicts import recursive_merge_dicts, sort_dict
from quacc.utils.files import get_uri
from quacc.wflow.db import results_to_db

if TYPE_CHECKING:
    from typing import Any

    from ase import Atoms
    from ase.calculators.calculator import Calculator
    from maggma.core import Store
    from phonopy import Phonopy

    from quacc.schemas._aliases.phonopy import PhononSchema


def summarize_phonopy(
    phonon: Phonopy,
    calculator: Calculator,
    input_atoms: Atoms | None = None,
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
    calculator
        Calculator used to generate the phonon object.
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
        "parameters": calculator.parameters,
        "phonopy_metadata": {"version": phonon.version},
        "nid": uri.split(":")[0],
        "dir_name": ":".join(uri.split(":")[1:]),
    }

    results = {"results": {"thermal_properties": phonon.get_thermal_properties_dict()}}
    phonon.save(settings={"force_constants": True})

    atoms_metadata = atoms_to_metadata(input_atoms) if input_atoms else {}
    unsorted_task_doc = recursive_merge_dicts(
        atoms_metadata, inputs, results, additional_fields
    )
    task_doc = sort_dict(unsorted_task_doc)

    if store:
        results_to_db(store, task_doc)

    return task_doc

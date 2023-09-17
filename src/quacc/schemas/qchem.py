"""Schemas for Q-Chem"""
from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

from maggma.core import Store
from monty.os.path import zpath
from pymatgen.io.qchem.inputs import QCInput
from pymatgen.io.qchem.outputs import QCOutput

from quacc import SETTINGS
from quacc.schemas.ase import summarize_run as base_summarize_run
from quacc.schemas.atoms import atoms_to_metadata
from quacc.utils.atoms import prep_next_run as prep_next_run_
from quacc.utils.db import results_to_db
from quacc.utils.dicts import clean_dict

if TYPE_CHECKING:
    from typing import TypeVar

    from ase import Atoms

    QchemSchema = TypeVar("QchemSchema")


def summarize_run(
    atoms: Atoms,
    dir_path: str | None = None,
    prep_next_run: bool = True,
    remove_empties: bool = False,
    additional_fields: dict | None = None,
    store: Store | None = None,
) -> QchemSchema:
    """
    Get tabulated results from a Q-Chem run and store them in a database-friendly
    format.

    Parameters
    ----------
    atoms
        ASE Atoms object following a calculation.
    dir_path
        Path to VASP outputs. A value of None specifies the current working
        directory
    prep_next_run
        Whether the Atoms object stored in {"atoms": atoms} should be prepared
        for the next run. This clears out any attached calculator and moves the
        final magmoms to the initial magmoms.
    remove_empties
        Whether to remove None values and empty lists/dicts from the
        TaskDocument.
    additional_fields
        Additional fields to add to the task document.
    store
        Maggma Store object to store the results in. If None,
        `SETTINGS.PRIMARY_STORE` will be used.

    Returns
    -------
    QchemSchema

        TODO: Description of fields.
    """

    additional_fields = additional_fields or {}
    dir_path = dir_path or Path.cwd()
    store = SETTINGS.PRIMARY_STORE if store is None else store

    base_summary = base_summarize_run(atoms, prep_next_run=prep_next_run, store=None)

    qc_output = {"qc_output": QCOutput(zpath(dir_path / "mol.qout")).data}
    qc_input = {"qc_input": QCInput.from_file(zpath(dir_path / "mol.qin")).as_dict()}

    task_doc = clean_dict(
        base_summary | qc_input | qc_output | additional_fields,
        remove_empties=remove_empties,
    )

    if store:
        results_to_db(store, task_doc)

    return task_doc

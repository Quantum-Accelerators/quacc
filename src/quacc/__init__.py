"""Init data for quacc."""

from __future__ import annotations

import logging
from importlib.metadata import version
from typing import TYPE_CHECKING

from ase.atoms import Atoms

from quacc.settings import QuaccSettings
from quacc.utils.dicts import Remove
from quacc.wflow_tools.customizers import redecorate, strip_decorator
from quacc.wflow_tools.decorators import Flow, Job, Subflow, flow, job, subflow

if TYPE_CHECKING:
    from typing import Any

__all__ = [
    "flow",
    "job",
    "subflow",
    "Flow",
    "Job",
    "Subflow",
    "redecorate",
    "strip_decorator",
    "Remove",
]


def atoms_as_dict(s: Atoms) -> dict[str, Any]:
    from ase.io.jsonio import encode
    from monty.json import jsanitize

    # Uses Monty's MSONable spec
    atoms_no_info = atoms.copy()
    atoms_no_info.info = {}
    return {"@module": "pymatgen.io.ase", "@class": "MSONAtoms", "atoms_json": encode(s),"atoms_no_info": jsanitize(atoms.info, strict=True)}


def atoms_from_dict(d: dict[str, Any]) -> Atoms:
    from ase.io.jsonio import decode
    from monty.json import MontyDecoder
    from pymatgen.io.ase import MSONAtoms

    # Uses Monty's MSONable spec
    mson_atoms = MSONAtoms(decode(dct["atoms_json"]))
    atoms_info = MontyDecoder().process_decoded(dct["atoms_info"])
    mson_atoms.info = atoms_info
    return mson_atoms


# Load the quacc version
__version__ = version("quacc")

# Make Atoms MSONable
Atoms.as_dict = atoms_as_dict
Atoms.from_dict = atoms_from_dict

# Load the settings
SETTINGS = QuaccSettings()

# Set logging info
logger = logging.getLogger(__name__)
level = logging.DEBUG if SETTINGS.DEBUG else logging.INFO
logger.setLevel(level)

# Monkeypatching for Prefect
if SETTINGS.WORKFLOW_ENGINE == "prefect":
    from prefect.futures import PrefectFuture
    from prefect.states import State

    def _patched_getitem(self, index):
        @job
        def _getitem(future, index_):
            return future[index_]

        return _getitem(self, index)

    PrefectFuture.__getitem__ = _patched_getitem
    State.__getitem__ = _patched_getitem

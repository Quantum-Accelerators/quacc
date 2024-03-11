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


def atoms_as_dict(atoms: Atoms) -> dict[str, Any]:
    """Uses Monty's MSONable spec. Refer to [pymatgen.io.ase.MSONAtoms][]"""

    from ase.io.jsonio import encode
    from monty.json import jsanitize

    atoms_no_info = atoms.copy()
    atoms_no_info.info = {}
    return {
        "@module": "pymatgen.io.ase",
        "@class": "MSONAtoms",
        "atoms_json": encode(atoms_no_info),
        "atoms_info": jsanitize(atoms.info, strict=True),
    }


def atoms_from_dict(dct: dict[str, Any]) -> Atoms:
    """Uses Monty's MSONable spec. Refer to [pymatgen.io.ase.MSONAtoms][]"""

    from ase.io.jsonio import decode
    from monty.json import MontyDecoder

    # Uses Monty's MSONable spec
    decoded_atoms = decode(dct["atoms_json"])
    atoms_info = MontyDecoder().process_decoded(dct["atoms_info"])
    decoded_atoms.info = atoms_info
    return decoded_atoms


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

    def _patched_getitem(self, index):
        @job
        def _getitem(future, index_):
            return future[index_]

        return _getitem(self, index)

    PrefectFuture.__getitem__ = _patched_getitem

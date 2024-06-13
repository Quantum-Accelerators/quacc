"""Init data for quacc."""

from __future__ import annotations

import logging
from importlib.metadata import version

from ase.atoms import Atoms
from pymatgen.io.ase import MSONAtoms

from quacc.settings import QuaccSettings, change_settings
from quacc.utils.dicts import Remove
from quacc.wflow_tools.customizers import redecorate, strip_decorator
from quacc.wflow_tools.decorators import Flow, Job, Subflow, flow, job, subflow

__all__ = [
    "flow",
    "job",
    "subflow",
    "Flow",
    "Job",
    "Subflow",
    "redecorate",
    "change_settings",
    "strip_decorator",
    "Remove",
]


# Load the quacc version
__version__ = version("quacc")

# Make Atoms MSONable
Atoms.as_dict = MSONAtoms.as_dict
Atoms.from_dict = MSONAtoms.from_dict

# Load the settings
SETTINGS = QuaccSettings()

# Set logging info
logging.basicConfig(level=logging.DEBUG if SETTINGS.DEBUG else logging.INFO)

# Monkeypatching for Prefect
if SETTINGS.WORKFLOW_ENGINE == "prefect":
    from prefect.client.schemas import State
    from prefect.futures import PrefectFuture

    def _patched_getitem(self, index):
        @job
        def _getitem(future, index_):
            return future[index_]

        return _getitem(self, index)

    PrefectFuture.__getitem__ = _patched_getitem
    State.__getitem__ = _patched_getitem

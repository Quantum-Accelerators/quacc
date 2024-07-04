"""Init data for quacc."""

from __future__ import annotations

import logging
import threading
from importlib.metadata import version
from typing import TYPE_CHECKING

from ase.atoms import Atoms
from pymatgen.io.ase import MSONAtoms

from quacc.settings import QuaccSettings, change_settings
from quacc.types import DefaultSetting
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
    "change_settings",
    "strip_decorator",
    "Remove",
    "get_settings",
    "QuaccDefault",
]


# Load the quacc version
__version__ = version("quacc")

# Make Atoms MSONable
Atoms.as_dict = MSONAtoms.as_dict  # type: ignore[attr-defined]
Atoms.from_dict = MSONAtoms.from_dict  # type: ignore[attr-defined]

# Load the settings
_thread_local = threading.local()


def _internally_set_settings(
    changes: dict[str, Any] | None = None, reset: bool = False
) -> None:
    """
    Set the `.settings` attribute for the current thread. This is not meant to be
    called by users. If you want to change the settings, use `from quacc import change_settings`

    Parameters
    ----------
    changes
        Changes to the settings, if any
    reset
        Reset the settings to the defaults of QuaccSettings()

    Returns
    -------
    None
    """
    changes = changes or {}
    if not hasattr(_thread_local, "settings") or reset:
        _thread_local.settings = QuaccSettings()
    for key, value in changes.items():
        setattr(_thread_local.settings, key, value)


def get_settings() -> QuaccSettings:
    """
    Get the `.settings` attribute for the current thread.

    Returns
    -------
    QuaccSettings
        The Pydantic settings model.
    """
    if not hasattr(_thread_local, "settings"):
        _internally_set_settings(reset=True)
    return _thread_local.settings


_settings = get_settings()

# Dummy value for when a default setting will be applied
QuaccDefault = DefaultSetting()

# Set logging info
logging.basicConfig(level=logging.DEBUG if _settings.DEBUG else logging.INFO)

# Monkeypatching for Prefect
if _settings.WORKFLOW_ENGINE == "prefect":
    from prefect.client.schemas import State
    from prefect.futures import PrefectFuture

    def _patched_getitem(self, index):
        @job
        def _getitem(future, index_):
            return future[index_]

        return _getitem(self, index)

    PrefectFuture.__getitem__ = _patched_getitem
    State.__getitem__ = _patched_getitem

"""Deprecated alias for quacc.recipes.mlip.phonons."""

from __future__ import annotations

import warnings

warnings.warn(
    "quacc.recipes.mlp.phonons is deprecated. Use quacc.recipes.mlip.phonons instead.",
    DeprecationWarning,
    stacklevel=2,
)

from quacc.recipes.mlip.phonons import *  # noqa: F403, E402
from quacc.recipes.mlip.phonons import phonon_flow  # noqa: F401, E402

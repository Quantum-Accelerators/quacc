"""Deprecated alias for quacc.recipes.mlip._base."""

import warnings

warnings.warn(
    "quacc.recipes.mlp._base is deprecated. Use quacc.recipes.mlip._base instead.",
    DeprecationWarning,
    stacklevel=2,
)

from quacc.recipes.mlip._base import *  # noqa: F401, F403, E402
from quacc.recipes.mlip._base import freezeargs, pick_calculator  # noqa: F401, E402

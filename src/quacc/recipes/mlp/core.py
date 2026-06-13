"""Deprecated alias for quacc.recipes.mlip.core."""

import warnings

warnings.warn(
    "quacc.recipes.mlp.core is deprecated. Use quacc.recipes.mlip.core instead.",
    DeprecationWarning,
    stacklevel=2,
)

from quacc.recipes.mlip.core import *  # noqa: F401, F403, E402
from quacc.recipes.mlip.core import relax_job, static_job  # noqa: F401, E402

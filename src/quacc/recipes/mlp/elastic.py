"""Deprecated alias for quacc.recipes.mlip.elastic."""

from __future__ import annotations

import warnings

warnings.warn(
    "quacc.recipes.mlp.elastic is deprecated. Use quacc.recipes.mlip.elastic instead.",
    DeprecationWarning,
    stacklevel=2,
)

from quacc.recipes.mlip.elastic import *  # noqa: F403, E402
from quacc.recipes.mlip.elastic import elastic_tensor_flow  # noqa: F401, E402

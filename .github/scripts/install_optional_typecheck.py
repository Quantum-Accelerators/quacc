"""Install optional API sources without GPU runtimes or model downloads.

Run after installing quacc[dev]. This checks signatures, not backend execution.
"""

from __future__ import annotations

import subprocess
import sys
import tomllib
from pathlib import Path

with Path("pyproject.toml").open("rb") as handle:
    extras = tomllib.load(handle)["project"]["optional-dependencies"]
requirements = sorted(
    {
        req
        for name, reqs in extras.items()
        if name not in {"dev", "docs"}
        for req in reqs
    }
)
# Extras are not resolved with --no-deps; these modules are imported directly.
requirements.extend(["distributed", "prefect-dask"])
subprocess.run(
    [sys.executable, "-m", "pip", "install", "--no-deps", *requirements], check=True
)
subprocess.run(
    [
        sys.executable,
        "-m",
        "pip",
        "install",
        "torch",
        "--index-url",
        "https://download.pytorch.org/whl/cpu",
    ],
    check=True,
)

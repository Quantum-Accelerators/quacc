from __future__ import annotations

from pathlib import Path

from quacc.recipes.espresso._base import prepare_copy
from quacc.wflow_tools.job_argument import Copy


def test_prepare_copy_preserves_copy_argument():
    copy_files = Copy({Path("previous-run"): ["charge-density.dat"]})

    assert prepare_copy(copy_files) is copy_files

from __future__ import annotations

import pytest

from quacc.utils.dependencies import requires


def test_requires_preserves_call_and_metadata():
    @requires(True, "missing")
    def calculate(value: int, *, scale: int = 2) -> int:
        """A guarded calculation."""
        return value * scale

    assert calculate(3, scale=4) == 12
    assert calculate.__name__ == "calculate"
    assert calculate.__doc__ == "A guarded calculation."


def test_requires_preserves_error():
    @requires(False, "Install the optional backend", ImportError)
    def calculate() -> None:
        pytest.fail("Missing dependencies must prevent execution")

    with pytest.raises(ImportError, match="Install the optional backend"):
        calculate()

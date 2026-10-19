"""Type-preserving guards for optional dependencies."""

from __future__ import annotations

from collections.abc import Callable
from typing import Any, TypeVar, cast

from monty.dev import requires as monty_requires

_F = TypeVar("_F", bound=Callable[..., Any])


def requires(
    condition: bool, message: str, err_cls: type[Exception] = RuntimeError
) -> Callable[[_F], _F]:
    """Guard a callable while preserving its signature.

    Delegates dependency checks and error messages to Monty. Its decorator
    currently returns a bare Callable, which otherwise loses parameter and
    return types on optional recipes.
    """
    guard = monty_requires(condition, message, err_cls)

    def decorate(function: _F) -> _F:
        return cast("_F", guard(function))

    return decorate

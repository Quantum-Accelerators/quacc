"""Execution context tracking for the autodiscover directory feature.

This module maintains a thread-/task-local stack of ContextNodes that records
the nesting of flows, subflows and jobs as they execute.  When
NESTED_RESULTS is enabled, the context stack is translated into a
filesystem path (e.g. ``my_flow-<id>/my_subflow-<id>/my_job-<id>``) so that
every job's output lands in a directory tree that mirrors the workflow
structure.

Two ContextVars are used:
  * ``_execution_context`` - the node stack (tuple of ContextNode).
  * ``_directory_context`` - the root directory chosen for the current
    top-level invocation.
"""

from __future__ import annotations

import asyncio
import shutil
from contextlib import contextmanager
from contextvars import ContextVar
from enum import Enum
from functools import wraps
from pathlib import Path
from typing import TYPE_CHECKING

from monty.json import MSONable

from quacc.utils.files import make_unique_name

if TYPE_CHECKING:
    from collections.abc import Callable


class NodeType(Enum):
    """Enumerates the types of nodes in a quacc workflow graph."""

    JOB = "job"
    FLOW = "flow"
    SUBFLOW = "subflow"


class ContextNode(MSONable):
    """A single entry in the execution-context stack.

    Each node records the function *name* and its *node_type* (job, flow, or
    subflow).  The stack of these nodes is used to build the autodiscovered
    directory path.  Inherits from MSONable so it can be serialized when the
    context is passed across process boundaries.
    """

    def __init__(self, name: str, node_type: NodeType):
        self.name = name
        self.node_type = node_type


# The node stack: tracks the current nesting of flow/subflow/job calls.
# An empty tuple means we are at the top level (no tracked function is active).
_execution_context: ContextVar[tuple[ContextNode, ...]] = ContextVar(
    "_execution_context", default=()
)

# The root output directory chosen for the current top-level invocation.
# Empty string means no directory has been assigned yet.
_directory_context: ContextVar[str] = ContextVar("_directory_context", default="")


def get_context() -> tuple[ContextNode, ...]:
    """Return the current execution-context node stack."""
    return _execution_context.get()


def get_directory_context() -> str:
    """Return the root output directory for the current top-level invocation."""
    return _directory_context.get()


def is_top_level() -> bool:
    """Check whether no tracked flow/job is currently executing."""
    current_context = get_context()
    return current_context == ()


def get_context_path() -> str:
    """Build a ``/``-separated path from the names of all nodes on the stack.

    For example, if the stack is ``[flow("bulk_to_slabs_flow"),
    subflow("bulk_to_slabs_subflow"), job("relax_job")]``, this returns
    ``"bulk_to_slabs_flow/bulk_to_slabs_subflow/relax_job"``.
    """
    return "/".join(n.name for n in get_context())


@contextmanager
def _push_context(name: str, node_type: NodeType):
    """Push a new node onto the execution-context stack for the duration of
    the ``with`` block, then restore the previous stack on exit."""
    old = _execution_context.get()
    new = (*old, ContextNode(name, node_type))
    token = _execution_context.set(new)
    try:
        yield
    finally:
        _execution_context.reset(token)


@contextmanager
def _push_directory_context(directory_path: str):
    """Temporarily set the root output directory for the current invocation."""
    token = _directory_context.set(directory_path)
    try:
        yield
    finally:
        _directory_context.reset(token)


def directory_context(directory_path: str):
    """Public convenience wrapper around ``_push_directory_context``."""
    return _push_directory_context(directory_path)


def tracked(node_type: NodeType):
    """
    Decorator factory that wraps a function to track its execution context.

    We use this to wrap functions before passing them to workflow engine decorators.

    For workflow engines that serialize tasks and execute them in separate
    processes (like Redun), the context can be passed via special kwargs
    (_quacc_ctx, _quacc_dir) which are extracted and restored before execution.
    """

    def decorator(func: Callable) -> Callable:
        # Preserve the true original (pre-decoration) function so that
        # strip_decorator() can fully unwrap back to the user's code.
        original = getattr(func, "original", func)

        if asyncio.iscoroutinefunction(func):

            @wraps(func)
            async def wrapper(*args, **kwargs):
                return await _tracked_call(func, node_type, args, kwargs)

        else:

            @wraps(func)
            def wrapper(*args, **kwargs):
                return _tracked_call(func, node_type, args, kwargs)

        # Attach the original unwrapped function so it can be recovered later.
        wrapper.original = original
        return wrapper

    return decorator


def _tracked_call(func, node_type, args, kwargs):
    """Shared tracking logic for both sync and async tracked wrappers.

    For async functions, the caller must ``await`` the return value.
    """
    from quacc import get_settings

    settings = get_settings()

    # When NESTED_RESULTS is off, skip all context tracking and just
    # delegate to the original function directly.
    if not settings.NESTED_RESULTS:
        return func(*args, **kwargs)

    # Create a unique name we can use at this level.
    name = make_unique_name(prefix=f"{func.__name__}-")

    if is_top_level():
        # This is the outermost tracked call: create a unique root
        # directory (e.g. ``<flow-name>-<timestamp>/``) and initialize both the
        # directory context and the execution-context stack.
        job_results_dir = settings.RESULTS_DIR.resolve()

        with directory_context(str(job_results_dir)), _push_context(name, node_type):
            return_value = func(*args, **kwargs)

            # If everything went okay, clean up temporary directory
            tmpdir_base = (settings.SCRATCH_DIR or settings.RESULTS_DIR).resolve()
            tmpdir = tmpdir_base / Path("tmp-" + name)
            if tmpdir.exists():
                shutil.rmtree(tmpdir)

            return return_value
    else:
        # We are inside an already-tracked invocation; just push another
        # node onto the existing stack.
        with _push_context(name, node_type):
            return func(*args, **kwargs)

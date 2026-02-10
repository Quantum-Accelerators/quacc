"""Workflow decorators."""

from __future__ import annotations

from collections.abc import Callable
from functools import partial, wraps
from typing import TYPE_CHECKING, Any

from quacc.settings import change_settings_wrap
from quacc.wflow_tools.context import (
    NodeType,
    get_context,
    get_directory_context,
    tracked,
)

if TYPE_CHECKING:
    from quacc.settings import QuaccSettings

Job = Callable[..., Any]
Flow = Callable[..., Any]
Subflow = Callable[..., Any]


def _make_context_capturing_wrapper(
    original_func: Callable, wrapped_callable: Callable
) -> Callable:
    """
    Create a wrapper that captures quacc context at call time and injects it
    for restoration at execution time.

    Parameters
    ----------
    original_func
        The original function being decorated (used for @wraps).
    wrapped_callable
        The workflow-engine-wrapped callable to invoke.

    Returns
    -------
    Callable
        A wrapper that captures context before calling the wrapped callable.
    """

    @wraps(original_func)
    def context_capturing_wrapper(*args, **kw):
        ctx = get_context()
        if ctx != ():
            kw["_quacc_ctx"] = ctx
        dir_ctx = get_directory_context()
        if dir_ctx != "":
            kw["_quacc_dir"] = dir_ctx
        return wrapped_callable(*args, **kw)

    return context_capturing_wrapper


def job(_func: Callable[..., Any] | None = None, **kwargs) -> Job:
    """
    Decorator for individual compute jobs. This is a `#!Python @job` decorator. Think of
    each `#!Python @job`-decorated function as an individual SLURM job, if that helps.

    | Quacc | Parsl        | Dask      | Prefect | Redun  | Jobflow |
    | ----- | ------------ | --------- | ------- | ------ | ------- |
    | `job` | `python_app` | `delayed` | `task`  | `task` | `job`   |

    All `#!Python @job`-decorated functions are transformed into their corresponding
    decorator.

    ```python
    from quacc import job


    @job
    def add(a, b):
        return a + b


    add(1, 2)
    ```

    ... is the same as doing

    === "Dask"

        ```python
        from dask import delayed


        @delayed
        def add(a, b):
            return a + b


        add(1, 2)
        ```

    === "Parsl"

        ```python
        from parsl import python_app


        @python_app
        def add(a, b):
            return a + b


        add(1, 2)
        ```

    === "Prefect"

        ```python
        from prefect import task


        @task
        def add(a, b):
            return a + b


        add.submit(1, 2)
        ```

    === "Redun"

        ```python
        from redun import task


        @task
        def add(a, b):
            return a + b


        add(1, 2)
        ```

    === "Jobflow"

        ```python
        import jobflow as jf


        @jf.job
        def add(a, b):
            return a + b


        add(1, 2)
        ```

    Parameters
    ----------
    _func
        The function to decorate. This is not meant to be supplied by the user.
    **kwargs
        Keyword arguments to pass to the workflow engine decorator.

    Returns
    -------
    Job
        The @job-decorated function.
    """
    from quacc import get_settings

    settings = get_settings()

    if _func is None:
        return partial(job, **kwargs)

    if changes := kwargs.pop("settings_swap", {}):
        _func = change_settings_wrap(_func, changes)

    if settings.WORKFLOW_ENGINE == "dask":
        from dask import delayed

        # Apply tracked() at runtime (inside the wrapper), not at decoration
        # time, because dask's tokenize/serialize can't handle the tracked
        # closure (it has a local __qualname__ and custom attributes that
        # trip up standard pickle).
        def wrapper(*f_args, **f_kwargs):
            return tracked(NodeType.JOB)(_func)(*f_args, **f_kwargs)

        wrapper.__name__ = _func.__name__
        wrapper.__wrapped__ = _func
        # Wrap in Delayed_ so the function is not eagerly evaluated when
        # composed inside a flow; Delayed_.__call__ captures context at
        # call time.
        return Delayed_(delayed(wrapper, **kwargs))

    # For non-Dask engines, wrap with tracked() now so every call pushes
    # a JOB node onto the execution-context stack.
    _func = tracked(NodeType.JOB)(_func)

    if settings.WORKFLOW_ENGINE == "jobflow":
        # Jobflow inspects function signatures to resolve references; add
        # explicit _quacc_ctx/_quacc_dir params so they survive inspection.
        _func_for_jf = _add_quacc_context_params(_func)
        _jf_func = _get_jobflow_wrapped_func(_func_for_jf, **kwargs)
        # Use the closure-based wrapper (not picklable, but jobflow doesn't
        # pickle the top-level callable).
        return _make_context_capturing_wrapper(_func, _jf_func)
    elif settings.WORKFLOW_ENGINE == "parsl":
        from parsl import python_app

        wrapped_fn = _get_parsl_wrapped_func(_func, kwargs)
        _parsl_func = python_app(wrapped_fn, **kwargs)
        # Use the class-based wrapper which is picklable (Parsl may serialize
        # decorated functions when passing them as arguments to other tasks).
        return _ContextCapturingWrapper(_func, _parsl_func)
    elif settings.WORKFLOW_ENGINE == "redun":
        from redun import task

        # Redun also inspects signatures, so add explicit context params.
        _func_for_redun = _add_quacc_context_params(_func)
        _redun_task = task(_func_for_redun, namespace=_func.__module__, **kwargs)
        # Class-based wrapper for pickle compatibility across Redun workers.
        return _ContextCapturingWrapper(_func, _redun_task)
    elif settings.WORKFLOW_ENGINE == "prefect":
        from prefect import task

        if settings.PREFECT_AUTO_SUBMIT:

            @wraps(_func)
            def wrapper(*f_args, **f_kwargs):
                decorated = task(_func, **kwargs)
                return decorated.submit(*f_args, **f_kwargs)

            return wrapper
        else:
            return task(_func, **kwargs)
    else:
        return _func


def flow(_func: Callable[..., Any] | None = None, **kwargs) -> Flow:
    """
    Decorator for workflows, which consist of at least one compute job. This is a
    `#!Python @flow` decorator.

    | Quacc  | Parsl     | Dask      | Prefect | Redun  | Jobflow   |
    | ------ | --------- | --------- | ------- | ------ | --------- |
    | `flow` | No effect | No effect | `flow`  | `task` | No effect |

    All `#!Python @flow`-decorated functions are transformed into their corresponding
    decorator.

    ```python
    from quacc import flow, job


    @job
    def add(a, b):
        return a + b


    @flow
    def workflow(a, b, c):
        return add(add(a, b), c)


    workflow(1, 2, 3)
    ```

    ... is the same as doing

    === "Dask"

        ```python
        from dask import delayed


        @delayed
        def add(a, b):
            return a + b


        def workflow(a, b, c):
            return add(add(a, b), c)


        workflow(1, 2, 3)
        ```

    === "Parsl"

        ```python
        from parsl import python_app


        @python_app
        def add(a, b):
            return a + b


        def workflow(a, b, c):
            return add(add(a, b), c)


        workflow(1, 2, 3)
        ```

    === "Prefect"

        ```python
        from prefect import flow, task


        @task
        def add(a, b):
            return a + b


        @flow
        def workflow(a, b, c):
            return add.submit(add.submit(a, b), c)


        workflow(1, 2, 3)
        ```

    === "Redun"

        ```python
        from redun import task


        @task
        def add(a, b):
            return a + b


        @task
        def workflow(a, b, c):
            return add(add(a, b), c)


        workflow(1, 2, 3)
        ```

    === "Jobflow"

        !!! Warning

            This decorator is not meant to be used with Jobflow at this time.

    Parameters
    ----------
    _func
        The function to decorate. This is not meant to be supplied by the user.
    **kwargs
        Keyword arguments to pass to the decorator.

    Returns
    -------
    Flow
        The `#!Python @flow`-decorated function.
    """
    from quacc import get_settings

    settings = get_settings()

    if _func is None:
        return partial(flow, **kwargs)

    # Wrap with tracked() so the flow pushes a FLOW node onto the
    # execution-context stack when it runs.
    _func = tracked(NodeType.FLOW)(_func)

    if settings.WORKFLOW_ENGINE == "redun":
        from redun import task

        _func_for_redun = _add_quacc_context_params(_func)
        _redun_task = task(_func_for_redun, namespace=_func.__module__, **kwargs)
        return _ContextCapturingWrapper(_func, _redun_task)
    elif settings.WORKFLOW_ENGINE == "prefect":
        return _get_prefect_wrapped_flow(_func, settings, **kwargs)
    elif settings.WORKFLOW_ENGINE == "jobflow":
        # Jobflow 0.3+ natively supports @flow, so delegate directly.
        return _get_jobflow_wrapped_flow(_func)
    elif settings.WORKFLOW_ENGINE == "parsl":
        # Parsl has no flow-level decorator; the tracked wrapper alone is
        # sufficient to establish context for child jobs.
        return _func
    else:
        return _func


def subflow(_func: Callable[..., Any] | None = None, **kwargs) -> Subflow:
    """
    Decorator for (dynamic) sub-workflows. This is a `#!Python @subflow` decorator.

    | Quacc     | Parsl      | Dask      | Prefect | Redun  | Jobflow   |
    | --------- | ---------- | --------- | ------- |------- | --------- |
    | `subflow` | `join_app` | `delayed` | `flow`  | `task` | No effect |

    All `#!Python @subflow`-decorated functions are transformed into their corresponding
    decorator.

    ```python
    import random
    from quacc import flow, job, subflow


    @job
    def add(a, b):
        return a + b


    @job
    def make_more(val):
        return [val] * random.randint(2, 5)


    @subflow
    def add_distributed(vals, c):
        return [add(val, c) for val in vals]


    @flow
    def workflow(a, b, c):
        result1 = add(a, b)
        result2 = make_more(result1)
        return add_distributed(result2, c)


    workflow(1, 2, 3)
    ```

    ... is the same as doing

    === "Dask"

        It's complicated... see the source code.

    === "Parsl"

        ```python
        import random
        from parsl import join_app, python_app


        @python_app
        def add(a, b):
            return a + b


        @python_app
        def make_more(val):
            return [val] * random.randint(2, 5)


        @join_app
        def add_distributed(vals, c):
            return [add(val, c) for val in vals]


        def workflow(a, b, c):
            result1 = add(a, b)
            result2 = make_more(result1)
            return add_distributed(result2, c)


        workflow(1, 2, 3)
        ```

    === "Prefect"

        ```python
        import random
        from prefect import flow, task


        @task
        def add(a, b):
            return a + b


        @task
        def make_more(val):
            return [val] * random.randint(2, 5)


        @flow
        def add_distributed(vals, c):
            return [add.submit(val, c) for val in vals]


        @flow
        def workflow(a, b, c):
            result1 = add.submit(a, b)
            result2 = make_more.submit(result1)
            return add_distributed(result2, c)


        workflow(1, 2, 3)
        ```

    === "Redun"

        ```python
        import random
        from redun import task


        @task
        def add(a, b):
            return a + b


        @task
        def make_more(val):
            return [val] * random.randint(2, 5)


        @task
        def add_distributed(vals, c):
            return [add(val, c) for val in vals]


        @task
        def workflow(a, b, c):
            result1 = add(a, b)
            result2 = make_more(result1)
            return add_distributed(result2, c)


        workflow(1, 2, 3)
        ```

    === "Jobflow"

        !!! Warning

            This decorator is not meant to be used with Jobflow at this time.

    Parameters
    ----------
    _func
        The function to decorate. This is not meant to be supplied by the user.
    **kwargs
        Keyword arguments to pass to the decorator.

    Returns
    -------
    callable
        The decorated function.
    """
    from quacc import get_settings

    settings = get_settings()

    if _func is None:
        return partial(subflow, **kwargs)

    if settings.WORKFLOW_ENGINE == "dask":
        from dask.delayed import delayed
        from dask.distributed import worker_client

        # Apply tracked() at runtime — same reason as @job above (dask
        # can't serialize the tracked closure at decoration time).
        def wrapper(*f_args, **f_kwargs):
            tracked_func = tracked(NodeType.SUBFLOW)(_func)
            # Use worker_client to compute nested delayed results inside
            # the already-running dask worker.
            with worker_client() as client:
                futures = client.compute(tracked_func(*f_args, **f_kwargs))
                return client.gather(futures)

        wrapper.__name__ = _func.__name__
        wrapper.__wrapped__ = _func
        return Delayed_(delayed(wrapper, **kwargs))

    # For non-Dask engines, wrap with tracked() now.
    _func = tracked(NodeType.SUBFLOW)(_func)

    if settings.WORKFLOW_ENGINE == "parsl":
        from parsl import join_app

        wrapped_fn = _get_parsl_wrapped_func(_func, kwargs)
        _parsl_func = join_app(wrapped_fn, **kwargs)
        return _ContextCapturingWrapper(_func, _parsl_func)
    elif settings.WORKFLOW_ENGINE == "prefect":
        return _get_prefect_wrapped_flow(_func, settings, **kwargs)
    elif settings.WORKFLOW_ENGINE == "redun":
        from redun import task

        _func_for_redun = _add_quacc_context_params(_func)
        _redun_task = task(_func_for_redun, namespace=_func.__module__, **kwargs)
        return _ContextCapturingWrapper(_func, _redun_task)
    elif settings.WORKFLOW_ENGINE == "jobflow":
        # Jobflow subflows are wrapped as @job (not @flow) so they appear
        # as replaceable jobs in the DAG.
        _func_for_jf = _add_quacc_context_params(_func)
        _jf_func = _get_jobflow_wrapped_func(_func_for_jf, **kwargs)
        return _make_context_capturing_wrapper(_func, _jf_func)
    else:
        return _func


class _ContextCapturingWrapper:
    """Picklable wrapper that captures quacc context at call time.

    Unlike the closure-based _make_context_capturing_wrapper, this class is
    picklable so that decorated functions can be passed as arguments to other
    tasks without serialization errors.  When pickled, it reduces to the
    underlying engine-wrapped callable.
    """

    def __init__(self, original_func, engine_callable):
        self._original_func = original_func
        self._engine_callable = engine_callable
        # Mirror standard function attributes so introspection works.
        self.__name__ = original_func.__name__
        self.__module__ = original_func.__module__
        self.__qualname__ = original_func.__qualname__
        self.__doc__ = original_func.__doc__
        self.__wrapped__ = original_func
        # Propagate .original from the tracked wrapper so that
        # strip_decorator() can unwrap all the way to the user's function.
        if hasattr(original_func, "original"):
            self.original = original_func.original

    def __call__(self, *args, **kw):
        # Snapshot the current context and inject it as kwargs so the
        # engine-wrapped callable (which may run in another process) can
        # restore the context on the worker side.
        ctx = get_context()
        if ctx != ():
            kw["_quacc_ctx"] = ctx
        dir_ctx = get_directory_context()
        if dir_ctx != "":
            kw["_quacc_dir"] = dir_ctx
        return self._engine_callable(*args, **kw)

    def __deepcopy__(self, memo):
        from copy import deepcopy

        # Preserve the wrapper (and context-capturing behavior) on deepcopy.
        # Functions aren't deepcopy'd (same object), but the engine callable
        # (e.g., redun Task) may be, so we deepcopy it.
        return _ContextCapturingWrapper(
            self._original_func, deepcopy(self._engine_callable, memo)
        )

    def __reduce__(self):
        # When serialized (e.g., passed as argument to another task),
        # serialize as the underlying engine callable so pickling succeeds.
        return self._engine_callable.__reduce__()


def _add_quacc_context_params(func: Callable) -> Callable:
    """Wrap *func* so ``_quacc_ctx`` and ``_quacc_dir`` appear as explicit
    keyword parameters in its signature.

    Jobflow and Redun both inspect function signatures to determine which
    arguments to resolve/pass.  By making the context params explicit
    (with defaults), the engines will accept them without error and
    forward them to the worker.
    """

    def wrapper(*args, _quacc_ctx=(), _quacc_dir="", **kwargs):
        # Re-inject context params into **kwargs so the tracked wrapper
        # (inside func) can pop them and restore the ContextVars.
        if _quacc_ctx != ():
            kwargs["_quacc_ctx"] = _quacc_ctx
        if _quacc_dir != "":
            kwargs["_quacc_dir"] = _quacc_dir
        return func(*args, **kwargs)

    # Copy function metadata so the engine sees the correct name/module.
    wrapper.__name__ = func.__name__
    wrapper.__module__ = func.__module__
    wrapper.__qualname__ = func.__qualname__
    wrapper.__doc__ = func.__doc__
    return wrapper


def _get_parsl_wrapped_func(
    func: Callable, decorator_kwargs: dict[str, Any]
) -> Callable:
    """
    Wrap a function to handle special Parsl arguments.

    Parameters
    ----------
    func
        The function to wrap.
    decorator_kwargs
        Decorator keyword arguments, including Parsl-specific ones that
        are meant to be passed to the underlying function. The `walltime`
        and `parsl_resource_specification` arguments will be injected
        into the function call and removed from the decorator arguments.

    Returns
    -------
    callable
        The wrapped function.
    """
    walltime = decorator_kwargs.pop("walltime", None)
    parsl_resource_specification = decorator_kwargs.pop(
        "parsl_resource_specification", {}
    )

    def wrapper(
        *f_args,
        walltime=walltime,  # noqa: ARG001
        parsl_resource_specification=parsl_resource_specification,  # noqa: ARG001
        # Context params added so Parsl's executor can forward them from the
        # _ContextCapturingWrapper on the submit side to the worker side.
        _quacc_ctx=(),
        _quacc_dir="",
        **f_kwargs,
    ):
        # Forward context into **kwargs so the tracked wrapper can restore it.
        if _quacc_ctx != ():
            f_kwargs["_quacc_ctx"] = _quacc_ctx
        if _quacc_dir != "":
            f_kwargs["_quacc_dir"] = _quacc_dir
        return func(*f_args, **f_kwargs)

    if getattr(func, "_changed", False):
        wrapper._changed = func._changed  # type: ignore[attr-defined]
        wrapper._original_func = func._original_func  # type: ignore[attr-defined]
    wrapper.__name__ = func.__name__
    return wrapper


def _get_prefect_wrapped_flow(
    _func: Callable, settings: QuaccSettings, **kwargs
) -> Callable:
    from prefect import flow as prefect_flow
    from prefect.futures import resolve_futures_to_results
    from prefect.utilities.asyncutils import is_async_fn

    if is_async_fn(_func):
        if settings.PREFECT_RESOLVE_FLOW_RESULTS:

            @wraps(_func)
            async def async_wrapper(*f_args, **f_kwargs):
                result = await _func(*f_args, **f_kwargs)
                return resolve_futures_to_results(result)

            return prefect_flow(async_wrapper, validate_parameters=False, **kwargs)

        else:
            return prefect_flow(_func, validate_parameters=False, **kwargs)
    else:
        if settings.PREFECT_RESOLVE_FLOW_RESULTS:

            @wraps(_func)
            def sync_wrapper(*f_args, **f_kwargs):
                result = _func(*f_args, **f_kwargs)
                return resolve_futures_to_results(result)

            return prefect_flow(sync_wrapper, validate_parameters=False, **kwargs)
        else:
            return prefect_flow(_func, validate_parameters=False, **kwargs)


def _get_jobflow_wrapped_func(method=None, **job_kwargs):
    from jobflow import job as jf_job

    return jf_job(method, **job_kwargs)


def _get_jobflow_wrapped_flow(_func: Callable) -> Callable:
    from jobflow import flow as jf_flow

    return jf_flow(_func)


class Delayed_:
    """A small Dask-compatible, serializable object to wrap delayed functions that we
    don't want to execute.

    At call time (in the flow's process), captures the current quacc context
    and injects it as kwargs so that the tracked wrapper on the worker can
    restore it.
    """

    __slots__ = ("func",)

    def __init__(self, func):
        self.func = func

    def __reduce__(self):
        return (Delayed_, (self.func,))

    def __call__(self, *args, **kwargs):
        # Capture the current quacc context at call time (in the flow's
        # process) and inject it as kwargs so that the tracked wrapper
        # running on the dask worker can restore the context stack.
        ctx = get_context()
        if ctx != ():
            kwargs["_quacc_ctx"] = ctx
        dir_ctx = get_directory_context()
        if dir_ctx != "":
            kwargs["_quacc_dir"] = dir_ctx
        return self.func(*args, **kwargs)

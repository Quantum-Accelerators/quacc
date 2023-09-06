from __future__ import annotations

import functools
from typing import TYPE_CHECKING

from monty.dev import requires

try:
    from prefect_dask.task_runners import DaskTaskRunner

    prefect_deps = True
except ImportError:
    prefect_deps = False
try:
    from dask_jobqueue import SLURMCluster

    dask_deps = True
except ImportError:
    dask_deps = False

if TYPE_CHECKING:
    from typing import Any, Callable, TypeVar

    from dask_jobqueue.core import Job as DaskJob

    Job = TypeVar("Job")
    Flow = TypeVar("Flow")
    Subflow = TypeVar("Subflow")


def job(_func: Callable | None = None, **kwargs) -> Job:  # sourcery skip
    """
    Decorator for individual compute jobs. This is a `@job` decorator. Think
    of each `@job`-decorated function as an individual SLURM job, if that helps.

    | Quacc | Covalent      | Parsl        | Jobflow | Redun  | Prefect |
    | ----- | ------------- | ------------ | ------- | ------ | ------- |
    | `job` | `ct.electron` | `python_app` | `job`   | `task` | `task`  |

    All `@job`-decorated functions are transformed into their corresponding
    decorator.

    ```python
    from quacc import job

    @job
    def add(a, b):
        return a + b
    ```

    is the same as doing

    ```python
    import covalent as ct

    @ct.electron
    def add(a, b):
        return a + b
    ```

    !!! Note

        For Prefect, this will actually call `Task.submit(*f_args, **f_kwargs)`, so you
        never need to call `.submit()` on your own.

    The wrapped function gets a new kwarg, `decorator_kwargs`, that can be used
    to modify the workflow engine decorator keyword arguments even after the
    quacc-decorated function has been imported. The wrapped (i.e. undecorated)
    function can also be stripped of its decorator by calling the `.__wrapped__`
    attribute.

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

    @functools.wraps(_func)
    def _inner(*f_args, decorator_kwargs: dict | None = None, **f_kwargs) -> Any:
        """
        The @job-decorated function.

        Parameters
        ----------
        *f_args
            Positional arguments to the function, if any.
        decorator_kwargs
            Keyword arguments to pass to the workflow engine decorator.
        **f_kwargs
            Keyword arguments to the function, if any.

        Returns
        -------
        Any
            The output of the @job-decorated function.
        """

        from quacc import SETTINGS

        if decorator_kwargs is None:
            decorator_kwargs = kwargs

        wflow_engine = SETTINGS.WORKFLOW_ENGINE
        if wflow_engine == "covalent":
            import covalent as ct

            decorated = ct.electron(_func, **decorator_kwargs)
        elif wflow_engine == "jobflow":
            from jobflow import job as jf_job

            decorated = jf_job(_func, **decorator_kwargs)
        elif wflow_engine == "parsl":
            from parsl import python_app

            decorated = python_app(_func, **decorator_kwargs)
        elif wflow_engine == "redun":
            from redun import task

            decorated = task(_func, **decorator_kwargs)
        elif wflow_engine == "prefect":
            from prefect import task

            decorated = task(_func, **decorator_kwargs)
        else:
            decorated = _func

        if wflow_engine == "prefect":
            return decorated.submit(*f_args, **f_kwargs)

        return decorated(*f_args, **f_kwargs)

    if _func is None:

        def decorator(_f):
            return job(_f, **kwargs)

        return decorator

    return _inner


def flow(_func: Callable | None = None, **kwargs) -> Flow:  # sourcery skip
    """
    Decorator for workflows, which consist of at least one compute job. This is
    a `@flow` decorator.

    | Quacc  | Covalent     | Parsl     | Jobflow   | Redun  | Prefect |
    | ------ | ------------ | --------- | --------- | ------ | ------- |
    | `flow` | `ct.lattice` | No effect | No effect | `task` | `flow`  |

    All `@flow`-decorated functions are transformed into their corresponding
    decorator.

    ```python
    from quacc import flow, job

    @job
    def add(a, b):
        return a + b

    @flow
    def workflow(a, b, c):
        return add(add(a, b), c)
    ```

    is the same as doing

    ```python
    import covalent as ct

    @ct.electron
    def add(a, b):
        return a + b

    @ct.lattice
    def workflow(a, b, c):
        return add(add(a, b), c)
    ```

    For the entries marked "No effect," the decorator has no effect on the
    underlying function.

    !!! Note

        For Covalent, this will actually call `ct.dispatch(ct.lattice(_func))(*f_args, **f_kwargs)`
        for the outermost `@flow` and just `ct.lattice(_func)(*f_args, **f_kwargs)`
        otherwise. As a result, you never need to call `ct.dispatch`.

    The wrapped function gets a new kwarg, `decorator_kwargs`, that can be used
    to modify the workflow engine decorator keyword arguments even after the
    quacc-decorated function has been imported. The wrapped (i.e. undecorated)
    function can also be stripped of its decorator by calling the `.__wrapped__`
    attribute.

    Parameters
    ----------
    _func
        The function to decorate. This is not meant to be supplied by the user.
    **kwargs
        Keyword arguments to pass to the decorator.

    Returns
    -------
    Flow
        The `@flow`-decorated function.
    """

    @functools.wraps(_func)
    def _inner(
        *f_args,
        decorator_kwargs: dict | None = None,
        dispatch_kwargs: dict | None = None,
        **f_kwargs,
    ) -> Any:
        """
        The @flow-decorated function.

        Parameters
        ----------
        *f_args
            Positional arguments to the function, if any.
        decorator_kwargs
            Keyword arguments to pass to the workflow engine decorator.
        dispatch_kwargs
            Keyword arguments to pass to `ct.dispatch()`, if Covalent is used.
        **f_kwargs
            Keyword arguments to the function, if any.

        Returns
        -------
        Any
            The output of the @flow-decorated function.
        """
        from quacc import SETTINGS

        if decorator_kwargs is None:
            decorator_kwargs = kwargs

        wflow_engine = SETTINGS.WORKFLOW_ENGINE

        dispatch_kwargs = dispatch_kwargs or {}
        if dispatch_kwargs and wflow_engine != "covalent":
            raise ValueError("The `dispatch_kwargs` argument only works with Covalent.")

        if wflow_engine == "covalent":
            import covalent as ct

            try:
                return ct.dispatch(ct.lattice(_func, **decorator_kwargs))(
                    *f_args, **f_kwargs
                )
            except (AttributeError, TypeError):
                return ct.lattice(_func, **dispatch_kwargs, **decorator_kwargs)(
                    *f_args, **f_kwargs
                )
        elif wflow_engine == "redun":
            from redun import task

            decorated = task(_func, **decorator_kwargs)
        elif wflow_engine == "prefect":
            from prefect import flow as prefect_flow

            if "validate_parameters" not in decorator_kwargs:
                decorator_kwargs["validate_parameters"] = False

            decorated = prefect_flow(_func, **decorator_kwargs)
        else:
            decorated = _func

        return decorated(*f_args, **f_kwargs)

    if _func is None:

        def decorator(_f):
            return flow(_f, **kwargs)

        return decorator

    return _inner


def subflow(_func: Callable | None = None, **kwargs) -> Subflow:  # sourcery skip
    """
    Decorator for (dynamic) sub-workflows. This is a `@subflow` decorator.

    | Quacc     | Covalent                  | Parsl      | Jobflow   | Redun  | Prefect |
    | --------- | ------------------------- | ---------- | --------- | ------ | ------- |
    | `subflow` | `ct.electron(ct.lattice)` | `join_app` | No effect | `task` | `flow`  |

    All `@subflow`-decorated functions are transformed into their corresponding
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
    ```

    is the same as doing

    ```python
    import covalent as ct
    import random

    @ct.electron
    def add(a, b):
        return a + b

    @ct.electron
    def make_more(val):
        return [val] * random.randint(2, 5)

    @ct.electron
    @ct.lattice
    def add_distributed(vals, c):
        return [add(val, c) for val in vals]

    @ct.lattice
    def workflow(a, b, c):
        result1 = add(a, b)
        result2 = make_more(result1)
        return add_distributed(result2, c)
    ```

    For the entries marked "No effect," the decorator has no effect on the
    underlying function.

    The wrapped function gets a new kwarg, `decorator_kwargs`, that can be used
    to modify the workflow engine decorator keyword arguments even after the
    quacc-decorated function has been imported. The wrapped (i.e. undecorated)
    function can also be stripped of its decorator by calling the `.__wrapped__`
    attribute.

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

    @functools.wraps(_func)
    def _inner(*f_args, decorator_kwargs: dict | None = None, **f_kwargs) -> Any:
        """
        The @subflow-decorated function.

        Parameters
        ----------
        *f_args
            Positional arguments to the function, if any.
        decorator_kwargs
            Keyword arguments to pass to the workflow engine decorator.
        **f_kwargs
            Keyword arguments to the function, if any.

        Returns
        -------
        Any
            The output of the @subflow-decorated function.
        """
        from quacc import SETTINGS

        if decorator_kwargs is None:
            decorator_kwargs = kwargs

        wflow_engine = SETTINGS.WORKFLOW_ENGINE
        if wflow_engine == "covalent":
            import covalent as ct

            decorated = ct.electron(ct.lattice(_func), **decorator_kwargs)
        elif wflow_engine == "parsl":
            from parsl import join_app

            decorated = join_app(_func, **decorator_kwargs)
        elif wflow_engine == "redun":
            from redun import task

            decorated = task(_func, **decorator_kwargs)
        elif wflow_engine == "prefect":
            from prefect import flow as prefect_flow

            if "validate_parameters" not in decorator_kwargs:
                decorator_kwargs["validate_parameters"] = False

            decorated = prefect_flow(_func, **decorator_kwargs)
        else:
            decorated = _func

        return decorated(*f_args, **f_kwargs)

    if _func is None:

        def decorator(_f):
            return subflow(_f, **kwargs)

        return decorator

    return _inner


@requires(prefect_deps and dask_deps, "Need quacc[prefect] dependencies")
def make_prefect_runner(
    cluster_kwargs: dict,
    cluster_class: callable = None,
    adapt_kwargs: dict[str, int | None] | None = None,
    client_kwargs: dict = None,
    temporary: bool = False,
) -> DaskTaskRunner:
    """
    Make a `DaskTaskRunner` for use with Prefect workflows.

    Parameters
    ----------
    cluster_kwargs
        Keyword arguments to pass to `cluster_class`.
    cluster_class
        The Dask cluster class to use. Defaults to `dask_jobqueue.SLURMCluster`.
    adapt_kwargs
        Keyword arguments to pass to `cluster.adapt` of the form `{"minimum": int, "maximum": int}`.
        If `None`, no adaptive scaling will be done.
    client_kwargs
        Keyword arguments to pass to `dask.distributed.Client`.
    temporary
        Whether to use a temporary cluster. If `True`, the cluster will be
        terminated once the `Flow` is finished. If `False`, the cluster will
        run until the walltime is reached and can run multiple `Flow`s.

    Returns
    -------
    DaskTaskRunner
        A DaskTaskRunner object for use with Prefect workflows.
    """

    if cluster_class is None:
        cluster_class = SLURMCluster

    # Make the one-time-use DaskTaskRunner
    if temporary:
        return DaskTaskRunner(
            cluster_class=cluster_class,
            cluster_kwargs=cluster_kwargs,
            adapt_kwargs=adapt_kwargs,
            client_kwargs=client_kwargs,
        )

    # Make the Dask cluster
    cluster = _make_dask_cluster(cluster_class, cluster_kwargs)

    # Set up adaptive scaling
    if adapt_kwargs and (adapt_kwargs["minimum"] or adapt_kwargs["maximum"]):
        cluster.adapt(minimum=adapt_kwargs["minimum"], maximum=adapt_kwargs["maximum"])

    # Return the DaskTaskRunner with the cluster address
    return DaskTaskRunner(address=cluster.scheduler_address)


@requires(dask_deps, "Need quacc[prefect] dependencies")
def _make_dask_cluster(
    cluster_class: callable, cluster_kwargs: dict, verbose=False
) -> DaskJob:
    """
    Make a Dask cluster for use with Prefect workflows.

    Parameters
    ----------
    cluster_class
        The Dask cluster class to use. Defaults to `dask_jobqueue.SLURMCluster`.
    cluster_kwargs
        Keyword arguments to pass to `cluster_class`.
    verbose
        Whether to print the job script to stdout.
    """
    cluster = cluster_class(**cluster_kwargs)
    if verbose:
        print(cluster.job_script())
    return cluster

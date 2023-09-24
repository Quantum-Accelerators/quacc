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


def job(_func: Callable | None = None, **kwargs) -> Job:
    """
    Decorator for individual compute jobs. This is a `#!Python @job` decorator. Think
    of each `#!Python @job`-decorated function as an individual SLURM job, if that helps.

    | Quacc | Covalent      | Parsl        | Prefect | Redun  | Jobflow |
    | ----- | ------------- | ------------ | ------- | ------ | ------- |
    | `job` | `ct.electron` | `python_app` | `task`  | `task` | `job`   |

    All `#!Python @job`-decorated functions are transformed into their corresponding
    decorator.

    The wrapped (i.e. undecorated) function can also be stripped of its decorator
    by calling the `#!Python .__wrapped__` attribute.

    ```python
    from quacc import job

    @job
    def add(a, b):
        return a + b

    add(1, 2)
    ```

    ... is the same as doing

    === "Covalent⭐"

        ```python
        import covalent as ct

        @ct.electron
        def add(a, b):
            return a + b

        add(1, 2)
        ```

    === "Parsl⭐"

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

    def decorator(func) -> Callable:
        @functools.wraps(func)
        def wrapper(f) -> Callable:
            @functools.wraps(f)
            def _inner(*f_args, **f_kwargs) -> Any:
                """
                This inner function is used for handling workflow engines that require some
                action beyond just decoration.

                Parameters
                ----------
                *f_args
                    Positional arguments to the function, if any.
                **f_kwargs
                    Keyword arguments to the function, if any.

                Returns
                -------
                Any
                    The output of the @job-decorated function.
                """
                if wflow_engine == "prefect":
                    return decorated_object.submit(*f_args, **f_kwargs)
                return decorated_object(*f_args, **f_kwargs)

            from quacc import SETTINGS

            wflow_engine = SETTINGS.WORKFLOW_ENGINE
            if wflow_engine == "covalent":
                import covalent as ct

                decorated_object = ct.electron(f, **kwargs)
            elif wflow_engine == "parsl":
                from parsl import python_app

                decorated_object = python_app(f, **kwargs)
            elif wflow_engine == "jobflow":
                import jobflow as jf

                decorated_object = jf.job(f, **kwargs)
            elif wflow_engine == "redun":
                from redun import task as redun_task

                decorated_object = redun_task(f, **kwargs)
            elif wflow_engine == "prefect":
                from prefect import task as prefect_task

                decorated_object = prefect_task(f, **kwargs)
                return _inner
            else:
                decorated_object = f

            if not hasattr(decorated_object, "__wrapped__"):
                decorated_object.__wrapped__ = _func

            return decorated_object

        return wrapper(func)

    if _func is None:
        return decorator
    return decorator(_func)


def flow(_func: Callable | None = None, **kwargs) -> Flow:
    """
    Decorator for workflows, which consist of at least one compute job. This is
    a `#!Python @flow` decorator.

    | Quacc  | Covalent     | Parsl     | Prefect | Redun  | Jobflow   |
    | ------ | ------------ | --------- | ------- | ------ | --------- |
    | `flow` | `ct.lattice` | No effect | `flow`  | `task` | No effect |

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

    === "Covalent⭐"

        ```python
        import covalent as ct

        @ct.electron
        def add(a, b):
            return a + b

        @ct.lattice
        def workflow(a, b, c):
            return add(add(a, b), c)

        workflow(1, 2, 3)
        ```

    === "Parsl⭐"

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

    def decorator(func) -> Callable:
        @functools.wraps(func)
        def wrapper(f) -> Callable:
            from quacc import SETTINGS

            wflow_engine = SETTINGS.WORKFLOW_ENGINE
            if wflow_engine == "covalent":
                import covalent as ct

                decorated_object = ct.lattice(f, **kwargs)
            elif wflow_engine == "redun":
                from redun import task as redun_task

                decorated_object = redun_task(f, **kwargs)
            elif wflow_engine == "prefect":
                from prefect import flow as prefect_flow

                decorated_object = prefect_flow(f, **kwargs)
            else:
                decorated_object = f
            return decorated_object

        return wrapper(func)

    if _func is None:
        return decorator
    return decorator(_func)


def subflow(_func: Callable | None = None, **kwargs) -> Subflow:
    """
    Decorator for (dynamic) sub-workflows. This is a `#!Python @subflow` decorator.

    | Quacc     | Covalent                  | Parsl      | Prefect | Redun  | Jobflow   |
    | --------- | ------------------------- | ---------- | ------- | ------ | --------- |
    | `subflow` | `ct.electron(ct.lattice)` | `join_app` | `flow`  | `task` | No effect |

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

    === "Covalent⭐"

        ```python
        import random
        import covalent as ct

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

        workflow(1, 2, 3)
        ```

    === "Parsl⭐"

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
            return [add(val, c) for val in vals]

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

    def decorator(func) -> Callable:
        @functools.wraps(func)
        def wrapper(f) -> Callable:
            from quacc import SETTINGS

            wflow_engine = SETTINGS.WORKFLOW_ENGINE
            if wflow_engine == "covalent":
                import covalent as ct

                decorated_object = ct.electron(ct.lattice(f, **kwargs))
            elif wflow_engine == "parsl":
                from parsl import join_app

                decorated_object = join_app(f, **kwargs)
            elif wflow_engine == "redun":
                from redun import task as redun_task

                decorated_object = redun_task(f, **kwargs)
            elif wflow_engine == "prefect":
                from prefect import flow as prefect_flow

                decorated_object = prefect_flow(f, **kwargs)
            else:
                decorated_object = f
            return decorated_object

        return wrapper(func)

    if _func is None:
        return decorator
    return decorator(_func)


@requires(prefect_deps and dask_deps, "Need quacc[prefect] dependencies")
def make_prefect_runner(
    cluster_kwargs: dict,
    cluster_class: Callable | None = None,
    adapt_kwargs: dict[str, int | None] | None = None,
    client_kwargs: dict | None = None,
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
    cluster_class: Callable, cluster_kwargs: dict, verbose=False
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
        print(cluster.job_script())  # noqa: T201
    return cluster

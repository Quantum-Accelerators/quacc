"""Workflow decorators."""
from __future__ import annotations

from functools import partial, wraps
from typing import TYPE_CHECKING, TypeVar

Job = TypeVar("Job")
Flow = TypeVar("Flow")
Subflow = TypeVar("Subflow")

if TYPE_CHECKING:
    from typing import Any, Callable


def job(_func: Callable | None = None, **kwargs) -> Job:
    """
    Decorator for individual compute jobs. This is a `#!Python @job` decorator. Think of
    each `#!Python @job`-decorated function as an individual SLURM job, if that helps.

    | Quacc | Covalent      | Parsl        | Dask      | Redun  | Jobflow |
    | ----- | ------------- | ------------ | --------- | ------ | ------- |
    | `job` | `ct.electron` | `python_app` | `delayed` | `task` | `job`   |

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

    === "Covalent"

        ```python
        import covalent as ct

        @ct.electron
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

    === "Dask"

        ```python
        from dask import delayed

        @delayed
        def add(a, b):
            return a + b

        add(1, 2)
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

    from quacc import SETTINGS

    if _func is None:
        return partial(job, **kwargs)

    elif SETTINGS.WORKFLOW_ENGINE == "covalent":
        import covalent as ct

        return ct.electron(_func, **kwargs)
    elif SETTINGS.WORKFLOW_ENGINE == "dask":
        from dask import delayed

        @wraps(_func)
        def wrapper(*args, **kwargs):
            return _func(*args, **kwargs)

        return Delayed_(delayed(wrapper))

    elif SETTINGS.WORKFLOW_ENGINE == "jobflow":
        from jobflow import job as jf_job

        return jf_job(_func, **kwargs)
    elif SETTINGS.WORKFLOW_ENGINE == "parsl":
        from parsl import python_app

        return python_app(_func, **kwargs)
    elif SETTINGS.WORKFLOW_ENGINE == "redun":
        from redun import task

        return task(_func, **kwargs)
    else:
        return _func


def flow(_func: Callable | None = None, **kwargs) -> Flow:
    """
    Decorator for workflows, which consist of at least one compute job. This is a
    `#!Python @flow` decorator.

    | Quacc  | Covalent     | Parsl     | Dask      | Redun  | Jobflow   |
    | ------ | ------------ | --------- | --------- | ------ | --------- |
    | `flow` | `ct.lattice` | No effect | No effect | `task` | No effect |

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

    === "Covalent"

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
    from quacc import SETTINGS

    if _func is None:
        return partial(flow, **kwargs)

    elif SETTINGS.WORKFLOW_ENGINE == "covalent":
        import covalent as ct

        return ct.lattice(_func, **kwargs)
    elif SETTINGS.WORKFLOW_ENGINE == "redun":
        from redun import task

        return task(_func, **kwargs)
    else:
        return _func


def subflow(_func: Callable | None = None, **kwargs) -> Subflow:
    """
    Decorator for (dynamic) sub-workflows. This is a `#!Python @subflow` decorator.

    | Quacc     | Covalent                  | Parsl      | Dask      | Redun  | Jobflow   |
    | --------- | ------------------------- | ---------- | --------- | ------ | --------- |
    | `subflow` | `ct.electron(ct.lattice)` | `join_app` | `delayed` | `task` | No effect |

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

    === "Covalent"

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

    === "Dask"

        ```python
        import random
        from dask import delayed

        @delayed
        def add(a, b):
            return a + b

        @delayed
        def make_more(val):
            return [val] * random.randint(2, 5)

        @delayed
        def add_distributed(vals, c):
            return [add(val, c) for val in vals]

        def workflow(a, b, c):
            result1 = add(a, b)
            result2 = make_more(result1)
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
    from quacc import SETTINGS

    if _func is None:
        return partial(subflow, **kwargs)

    elif SETTINGS.WORKFLOW_ENGINE == "covalent":
        import covalent as ct

        return ct.electron(ct.lattice(_func), **kwargs)
    elif SETTINGS.WORKFLOW_ENGINE == "parsl":
        from parsl import join_app

        return join_app(_func, **kwargs)
    elif SETTINGS.WORKFLOW_ENGINE == "redun":
        from redun import task

        return task(_func, **kwargs)
    elif SETTINGS.WORKFLOW_ENGINE == "dask":
        from dask import delayed
        from dask.distributed import worker_client

        @wraps(_func)
        def wrapper(*args, **kwargs):
            with worker_client() as client:
                futures = client.compute(_func(*args, **kwargs))
                results = client.gather(futures)
                return results

        return delayed(wrapper)
    else:
        return _func


class Delayed_:
    """
    A small Dask-compatible, serializable object to wrap delayed functions
    that we don't want to execute
    """

    __slots__ = ("func",)

    def __init__(self, func):
        self.func = func

    def __reduce__(self):
        return (Delayed_, (self.func,))

    def __call__(self, *args, **kwargs):
        return self.func(*args, **kwargs)

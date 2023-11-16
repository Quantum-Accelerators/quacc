"""Workflow decorators."""
from __future__ import annotations

import functools
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

    | Quacc | Covalent      | Parsl        | Prefect | Redun  | Jobflow |
    | ----- | ------------- | ------------ | ------- | ------ | ------- |
    | `job` | `ct.electron` | `python_app` | `task`  | `task` | `job`   |

    All `#!Python @job`-decorated functions are transformed into their corresponding
    decorator.

    The wrapped function can also be stripped of its decorator by calling the
    `#!Python .__wrapped__` attribute.

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

    @functools.wraps(_func)
    def _inner(
        *f_args, decorator_kwargs: dict[str, Any] | None = None, **f_kwargs
    ) -> Any:
        """
        This function is used for handling workflow engines that require some action
        beyond just decoration. It also patches the parent function `_func` to takke an
        additional keyword argument, `deocrator_kwargs`, that is a dictionary of keyword
        arguments to pass during the decorator construction.

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
        decorator_kwargs = decorator_kwargs if decorator_kwargs is not None else kwargs

        if wflow_engine == "prefect":
            from prefect import task as prefect_task

            decorated = prefect_task(_func, **decorator_kwargs)
            return decorated.submit(*f_args, **f_kwargs)

        return decorated(*f_args, **f_kwargs)

    from quacc import SETTINGS

    wflow_engine = SETTINGS.WORKFLOW_ENGINE

    if _func is None:
        return functools.partial(job, **kwargs)

    if wflow_engine == "covalent":
        import covalent as ct

        decorated = ct.electron(_func, **kwargs)
    elif wflow_engine == "jobflow":
        from jobflow import job as jf_job

        decorated = jf_job(_func, **kwargs)
    elif wflow_engine == "parsl":
        from parsl import python_app

        decorated = python_app(_func, **kwargs)
    elif wflow_engine == "redun":
        from redun import task as redun_task

        decorated = redun_task(_func, **kwargs)
    elif wflow_engine == "prefect":
        return _inner
    else:
        decorated = _func

    if not hasattr(decorated, "__wrapped__"):
        decorated.__wrapped__ = _func

    return decorated


def flow(
    _func: Callable | None = None, **kwargs
) -> Flow:  # sourcery skip: lift-return-into-if, switch
    """
    Decorator for workflows, which consist of at least one compute job. This is a
    `#!Python @flow` decorator.

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
    from quacc import SETTINGS

    if _func is None:
        return functools.partial(flow, **kwargs)

    wflow_engine = SETTINGS.WORKFLOW_ENGINE
    if wflow_engine == "covalent":
        import covalent as ct

        decorated = ct.lattice(_func, **kwargs)
    elif wflow_engine == "prefect":
        from prefect import flow as prefect_flow

        decorated = prefect_flow(_func, **kwargs)
    elif wflow_engine == "redun":
        from redun import task as redun_task

        decorated = redun_task(_func, **kwargs)
    else:
        decorated = _func

    return decorated


def subflow(
    _func: Callable | None = None, **kwargs
) -> Subflow:  # sourcery skip: lift-return-into-if, switch
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

    from quacc import SETTINGS

    if _func is None:
        return functools.partial(subflow, **kwargs)

    wflow_engine = SETTINGS.WORKFLOW_ENGINE
    if wflow_engine == "covalent":
        import covalent as ct

        decorated = ct.electron(ct.lattice(_func, **kwargs))
    elif wflow_engine == "parsl":
        from parsl import join_app

        decorated = join_app(_func, **kwargs)
    elif wflow_engine == "prefect":
        from prefect import flow as prefect_flow

        decorated = prefect_flow(_func, **kwargs)
    elif wflow_engine == "redun":
        from redun import task as redun_task

        decorated = redun_task(_func, **kwargs)
    else:
        decorated = _func

    return decorated

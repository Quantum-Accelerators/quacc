from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from typing import TypeVar

    from ase import Atoms

    Job = TypeVar("Job")
    Flow = TypeVar("Flow")
    Subflow = TypeVar("Subflow")


def job(_func: callable | None = None, **kwargs) -> Job:  # sourcery skip
    """
    Decorator for individual compute jobs. This is a @job decorator.

    @job = @ct.electron [Covalent] = @python_app [Parsl] = @job [Jobflow]

    Parameters
    ----------
    _func
        The function to decorate.
    **kwargs
        Keyword arguments to pass to the decorator.

    Returns
    -------
    Job
        The @job-decorated function. The decorated function will have an attribute `__wrapped__`
        which is the original (unwrapped) function.
    """

    from quacc import SETTINGS

    wflow_engine = (
        SETTINGS.WORKFLOW_ENGINE.lower() if SETTINGS.WORKFLOW_ENGINE else None
    )
    if wflow_engine == "covalent":
        import covalent as ct

        decorated = ct.electron(_func, **kwargs)
    elif wflow_engine == "jobflow":
        from jobflow import job as jf_job

        decorated = jf_job(_func, **kwargs)
    elif wflow_engine == "parsl":
        from parsl import python_app

        decorated = python_app(_func, **kwargs)
    else:
        decorated = _func

    decorated.__wrapped__ = _func

    return decorated


def flow(_func: callable | None = None, **kwargs) -> Flow:  # sourcery skip
    """
    Decorator for workflows, which consist of at least one compute job. This is a @flow decorator.

    @flow = @ct.lattice [Covalent]. For Parsl and Jobflow, the decorator returns the
    original function, unchanged. The decorated function will have an attribute `__wrapped__`
    which is the original (unwrapped) function.

    Parameters
    ----------
    _func
        The function to decorate.
    **kwargs
        Keyword arguments to pass to the decorator.

    Returns
    -------
    Flow
        The @flow-decorated function.
    """

    from quacc import SETTINGS

    wflow_engine = (
        SETTINGS.WORKFLOW_ENGINE.lower() if SETTINGS.WORKFLOW_ENGINE else None
    )
    if wflow_engine == "covalent":
        import covalent as ct

        decorated = ct.lattice(_func, **kwargs)
    else:
        decorated = _func

    return decorated


def subflow(_func: callable | None = None, **kwargs) -> Subflow:  # sourcery skip
    """
    Decorator for (dynamic) sub-workflows. This is a @subflow decorator.

    @subflow = @ct.electron(@ct.lattice) [Covalent] = @join_app [Parsl].
    For Jobflow, the decorator returns the original (unwrapped) function.
    The decorated function will have an attribute `__wrapped__`
    which is the original (unwrapped) function.

    Parameters
    ----------
    _func
        The function to decorate.
    **kwargs
        Keyword arguments to pass to the decorator.

    Returns
    -------
    callable
        The decorated function.
    """

    from quacc import SETTINGS

    wflow_engine = (
        SETTINGS.WORKFLOW_ENGINE.lower() if SETTINGS.WORKFLOW_ENGINE else None
    )
    if wflow_engine == "covalent":
        import covalent as ct

        decorated = ct.electron(ct.lattice(_func), **kwargs)
    elif wflow_engine == "parsl":
        from parsl import join_app

        decorated = join_app(_func, **kwargs)
    else:
        decorated = _func

    return decorated


def fetch_atoms(atoms: Atoms | dict) -> Atoms:
    """
    Returns an Atoms object from a typical quacc recipe input, which can
    either be an `Atoms` object or a dictionary with an entry `{"atoms": Atoms}`.
    It may seem a bit silly to do this, but there is a purpose. If you want to
    write a workflow where the output of one recipe is passed to the input of
    another recipe, you can always do output["atoms"] to fetch the output Atoms
    object to pass to the input to the second function. However, this process
    will often be its own compute step in workflow management tools because they
    need to resolve the output in order to query it. Depending on the workflow manager,
    this can be a waste of compute resources, so it's oftentimes better to do this
    parsing inside the compute task itself, which is why passing in the raw dictionary
    can be preferred.

    Parameters
    ----------
    atoms
        Atoms object or dictionary with an entry {"atoms": Atoms}

    Returns
    -------
    Atoms
        Atoms object
    """
    try:
        return atoms["atoms"]
    except Exception:
        return atoms

from typing import Any
from ase.atoms import Atom, Atoms
from ase.io.jsonio import encode, decode
from monty.json import jsanitize as jsanitize_
from monty.json import MSONable, MontyDecoder
import numpy as np


def jsanitize(obj: Any) -> Any:
    """
    Sanitize an input object such that it is JSON serializable. This is
    a wrapper around Monty's jsanitize function but is modified to deal
    with Atom/Atoms objects and better handling of MSONables.

    Anything being passed as an input or output for Jobflow should be
    jsanitized.

    Parameters
    ----------
    obj
        Object to sanitize.

    Returns
    -------
    Any
        The JSON-serializable object.
    """

    if isinstance(obj, (Atom, Atoms)):
        return encode(obj)
    if isinstance(obj, (list, tuple, np.ndarray)):
        return [jsanitize(o) for o in obj]
    if isinstance(obj, dict):
        return {k.__str__(): jsanitize(v) for k, v in obj.items()}
    if isinstance(obj, MSONable):
        return obj.as_dict()

    return jsanitize_(obj)


def jdesanitize(obj: Any) -> Any:
    """
    Desanitize an input object that was jsanitized. This is
    a wrapper around Monty's jsanitize function but is modified to deal
    with Atom/Atoms objects and better handling of MSONables.

    Anything being passed as an input or output for Jobflow should be
    jsanitized.

    Parameters
    ----------
    obj
        Object to sanitize.

    Returns
    -------
    Any
        The JSON-serializable object.
    """
    if isinstance(obj, str) and (
        '"__ase_objtype__": "atoms"' in obj or '"__ase_objtype__": "atom"' in obj
    ):
        return decode(obj)
    if isinstance(obj, list):
        return [jdesanitize(o) for o in obj]
    if isinstance(obj, dict):
        if "@module" in obj.keys():
            return MontyDecoder().process_decoded(obj)
        return {k.__str__(): jdesanitize(v) for k, v in obj.items()}

    return obj

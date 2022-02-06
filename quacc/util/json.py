from typing import Any

import numpy as np
from ase.atoms import Atom, Atoms
from ase.io.jsonio import decode, encode
from monty.json import MontyDecoder, MSONable, jsanitize


def jsonify(obj: Any) -> Any:
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
        return [jsonify(o) for o in obj]
    if isinstance(obj, dict):
        return {k.__str__(): jsonify(v) for k, v in obj.items()}
    if isinstance(obj, MSONable):
        return obj.as_dict()

    return jsanitize(obj)


def unjsonify(obj: Any) -> Any:
    """
    Desanitize an input object that was jsanitized. This is
    a wrapper around Monty's MontyDecoer() function but is modified to deal
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
        return [unjsonify(o) for o in obj]
    if isinstance(obj, dict):
        if "@module" in obj.keys():
            return MontyDecoder().process_decoded(obj)
        return {k.__str__(): unjsonify(v) for k, v in obj.items()}

    return obj

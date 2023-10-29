"""
Parameter-related utilities for the Q-Chem calculator.
"""
from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from typing import Any


def cleanup_params(
    qchem_input_params: dict[str, Any], method: str | None
) -> dict[str, Any]:
    """
    Cleanup parameters for the Q-Chem calculator.

    Parameters
    ----------
    qchem_input_params
        Dictionary of Q-Chem input parameters to be passed to
        `pymatgen.io.qchem.sets.DictSet`.
    method
        The method to use for the calculation.

    Returns
    -------
    dict[str, Any]
        The cleaned up Q-Chem input parameters.
    """

    if "overwrite_inputs" not in qchem_input_params:
        qchem_input_params["overwrite_inputs"] = {}

    if qchem_input_params.get("smd_solvent") and qchem_input_params.get(
        "pcm_dielectric"
    ):
        raise ValueError("PCM and SMD cannot be employed simultaneously.")

    if "rem" not in qchem_input_params["overwrite_inputs"]:
        qchem_input_params["overwrite_inputs"]["rem"] = {}
    if method and "method" not in qchem_input_params["overwrite_inputs"]["rem"]:
        qchem_input_params["overwrite_inputs"]["rem"]["method"] = method
    return qchem_input_params


def get_default_params(
    cores: int,
    charge: int,
    spin_multiplicity: int,
    method: str | None,
    basis_set: str,
    scf_algorithm: str,
    qchem_input_params: dict[str, Any],
) -> dict[str, Any]:
    """
    Get the parameters that have been passed to the Q-Chem
    calculator to store later in FileIOCalculator's self.default_parameters

    Parameters
    ----------
    cores
        Number of cores to use for the Q-Chem calculation.
    charge
        The total charge of the molecular system.
    spin_multiplicity
        The spin multiplicity of the molecular system.
    method
        The method to use for the calculation.
    basis_set
        The basis set to use for the calculation.
    scf_algorithm
        The SCF algorithm to use for the calculation.
    qchem_input_params
        Dictionary of Q-Chem input parameters to be passed to
        `pymatgen.io.qchem.sets.DictSet`.

    Returns
    -------
    dict[str, Any]
        The default parameters that have been passed to the Q-Chem calculator.
    """
    default_parameters = {
        "cores": cores,
        "charge": charge,
        "spin_multiplicity": spin_multiplicity,
        "scf_algorithm": scf_algorithm,
        "basis_set": basis_set,
    }

    if method:
        default_parameters["method"] = method

    # We also want to save the contents of self.qchem_input_params. However,
    # the overwrite_inputs key will have a corresponding value which is
    # either an empty dictionary or a nested dict of dicts, requiring a bit
    # of careful unwrapping.
    for key in qchem_input_params:
        if key == "overwrite_inputs":
            for subkey in qchem_input_params[key]:
                for subsubkey in qchem_input_params[key][subkey]:
                    default_parameters[
                        f"overwrite_{subkey}_{subsubkey}"
                    ] = qchem_input_params[key][subkey][subsubkey]
        else:
            default_parameters[key] = qchem_input_params[key]

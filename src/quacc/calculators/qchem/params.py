"""Parameter-related utilities for the Q-Chem calculator."""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING

from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.io.qchem.inputs import QCInput
from pymatgen.io.qchem.sets import QChemDictSet
from pymatgen.io.qchem.utils import lower_and_check_unique

from quacc.utils.dicts import recursive_dict_merge, sort_dict

if TYPE_CHECKING:
    from typing import Any

    from ase.atoms import Atoms

    from quacc.calculators.qchem.qchem import QChem

logger = logging.getLogger(__name__)


def make_qc_input(qchem: QChem, atoms: Atoms) -> QCInput:
    """
    Make a QCInput object. It will, by default, create a QCInput from the QChem
    calculator kwargs. If `qchem.qchem_dict_set_params` is specified, it will create a
    QCInput from a QChemDictSet, merging the two QCInput objects and taking the latter
    as higher priority.

    Parameters
    ----------
    qchem
        The QChem object.

    Returns
    -------
    QCInput
        The QCInput object.
    """
    atoms.charge = qchem.charge  # type: ignore[attr-defined]
    atoms.spin_multiplicity = qchem.spin_multiplicity  # type: ignore[attr-defined]
    molecule = AseAtomsAdaptor().get_molecule(atoms)

    if qchem.qchem_dict_set_params:
        # Get minimal parameters needed to instantiate a QChemDictSet
        if "molecule" in qchem.qchem_dict_set_params:
            msg = "Do not specify `molecule` in `qchem_dict_set_params`"
            raise NotImplementedError(msg)
        if "job_type" not in qchem.qchem_dict_set_params and qchem.rem.get("job_type"):
            qchem.qchem_dict_set_params["job_type"] = qchem.rem["job_type"]
        if "basis_set" not in qchem.qchem_dict_set_params and qchem.rem.get("basis"):
            qchem.qchem_dict_set_params["basis_set"] = qchem.rem["basis"]
        if "scf_algorithm" not in qchem.qchem_dict_set_params and qchem.rem.get(
            "scf_algorithm"
        ):
            qchem.qchem_dict_set_params["scf_algorithm"] = qchem.rem["scf_algorithm"]
        if "qchem_version" not in qchem.qchem_dict_set_params:
            qchem.qchem_dict_set_params["qchem_version"] = 6

        # Make QChemDictSet
        qc_dict_set = QChemDictSet(molecule, **qchem.qchem_dict_set_params)
        for prop in [
            "rem",
            "opt",
            "pcm",
            "solvent",
            "smx",
            "scan",
            "van_der_waals",
            "plots",
            "nbo",
            "geom_opt",
            "svp",
            "pcm_nonels",
        ]:
            prop1 = getattr(qchem, prop)
            if prop2 := getattr(qc_dict_set, prop):
                setattr(qchem, prop, recursive_dict_merge(prop2, prop1))
        for prop in ["vdw_mode", "cdft", "almo_coupling"]:
            prop2 = getattr(qc_dict_set, prop)
            if prop2 and not prop1:
                setattr(qchem, prop, prop2)

    qchem.rem = sort_dict(
        get_rem_swaps(qchem.rem, restart=qchem.prev_orbital_coeffs is not None)
    )

    return QCInput(
        molecule,
        qchem.rem,
        opt=qchem.opt,
        pcm=qchem.pcm,
        solvent=qchem.solvent,
        smx=qchem.smx,
        scan=qchem.scan,
        van_der_waals=qchem.van_der_waals,
        vdw_mode=qchem.vdw_mode,
        plots=qchem.plots,
        nbo=qchem.nbo,
        geom_opt=qchem.geom_opt,
        cdft=qchem.cdft,
        almo_coupling=qchem.almo_coupling,
        svp=qchem.svp,
        pcm_nonels=qchem.pcm_nonels,
    )


def cleanup_attrs(qchem: QChem) -> None:
    """
    Clean up self attribute parameters in place.

    Parameters
    ----------
    qchem
        The QChem object

    Returns
    -------
    None
    """
    for attr in [
        "rem",
        "pcm",
        "solvent",
        "smx",
        "scan",
        "van_der_waals",
        "plots",
        "nbo",
        "geom_opt",
        "svp",
        "pcm_nonels",
    ]:
        attr_val = lower_and_check_unique(getattr(qchem, attr))
        setattr(qchem, attr, attr_val)


def get_rem_swaps(rem: dict[str, Any], restart: bool = False) -> dict[str, Any]:
    """
    Automatic swaps for the rem dictionary.

    Parameters
    ----------
    rem
        rem dictionary
    restart
        Whether or not this is a restart calculation.

    Returns
    -------
    dict[str, Any]
        rem dictionary with swaps
    """
    if restart and "scf_guess" not in rem:
        logger.info("Copilot: Setting scf_guess in `rem` to 'read'")
        rem["scf_guess"] = "read"
    if "max_scf_cycles" not in rem:
        rem["max_scf_cycles"] = 200 if rem.get("scf_algorithm") == "gdm" else 100
        logger.info(
            f"Copilot: Setting max_scf_cycles in `rem` to {rem['max_scf_cycles']}"
        )

    return rem

"""Q-Chem calculator IO"""
from __future__ import annotations

import struct
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np
from ase import units
from ase.calculators.calculator import FileIOCalculator
from emmet.core.tasks import _parse_custodian
from monty.io import zopen
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.io.qchem.inputs import QCInput
from pymatgen.io.qchem.outputs import QCOutput
from pymatgen.io.qchem.sets import QChemDictSet

if TYPE_CHECKING:
    from typing import Any, Literal, TypedDict

    from ase import Atoms
    from numpy.typing import NDArray

    class Results(TypedDict, total=False):
        energy: float  # electronic energy in eV
        forces: NDArray  # forces in eV/A
        hessian: NDArray  # Hessian in native Q-Chem units
        enthalpy: float  # total enthalpy in eV
        entropy: float  # total entropy in eV/K
        qc_output: dict[
            str, Any
        ]  # Output from `pymatgen.io.qchem.outputs.QCOutput.data`
        qc_input: dict[
            str, Any
        ]  # Input from `pymatgen.io.qchem.inputs.QCInput.as_dict()`
        custodian: dict[str, Any]  # custodian.json file metadata


def write_qchem(
    atoms: Atoms,
    charge: int = 0,
    spin_multiplicity: int = 1,
    basis_set: str = "def2-tzvpd",
    job_type: Literal = ["sp", "force", "opt", "freq"],
    scf_algorithm: str = "diis",
    qchem_input_params: dict[str, Any] | None = None,
    prev_orbital_coeffs: list[float] | None = None,
    properties: list[str] | None = None,
    system_changes: list[str] | None = None,
) -> None:
    """
    Write the Q-Chem input files.

    Parameters
    ----------
    atoms
        The Atoms object to be used for the calculation.
    charge
        The total charge of the molecular system.
    spin_multiplicity
        The spin multiplicity of the molecular system.
    basis_set
        The basis set to use for the calculation.
    job_type
        The type of calculation to perform.
    scf_algorithm
        The SCF algorithm to use for the calculation.
    qchem_input_params
        Dictionary of Q-Chem input parameters to be passed to
        `pymatgen.io.qchem.sets.DictSet`.
    prev_orbital_coeffs
        The orbital coefficients from a previous calculation.
    properties
        List of properties to calculate.
    system_changes
        List of system changes to make.

    Returns
    -------
    None
    """
    FileIOCalculator.write_input(atoms, properties, system_changes)

    atoms.charge = charge
    atoms.spin_multiplicity = spin_multiplicity
    mol = AseAtomsAdaptor.get_molecule(atoms)

    if prev_orbital_coeffs is not None:
        with Path("53.0").open(mode="wb") as file:
            for val in prev_orbital_coeffs:
                data = struct.pack("d", val)
                file.write(data)
        if "overwrite_inputs" not in qchem_input_params:
            qchem_input_params["overwrite_inputs"] = {}
        if "rem" not in qchem_input_params["overwrite_inputs"]:
            qchem_input_params["overwrite_inputs"]["rem"] = {}
        if "scf_guess" not in qchem_input_params["overwrite_inputs"]["rem"]:
            qchem_input_params["overwrite_inputs"]["rem"]["scf_guess"] = "read"

    qcin = QChemDictSet(
        mol,
        job_type,
        basis_set,
        scf_algorithm,
        qchem_version=6,
        **qchem_input_params,
    )
    qcin.write("mol.qin")


def read_qchem(
    job_type: Literal = ["sp", "force", "opt", "freq"]
) -> tuple[Results, list[float]]:
    """
    Read Q-Chem log files.

    Parameters
    ----------
    job_type
        The type of calculation to perform.

    Returns
    -------
    tuple[Results, list[float]]
        The results of the calculation and the orbital coefficients from a previous
        calculation.
    """
    qc_input = QCInput.from_file("mol.qin").as_dict()
    qc_output = QCOutput("mol.qout").data

    results: Results = {}

    results["energy"] = qc_output["final_energy"] * units.Hartree
    results["qc_output"] = qc_output
    results["qc_input"] = qc_input
    results["custodian"] = _parse_custodian(Path.cwd())

    if job_type in ["force", "opt"]:
        # Read the gradient scratch file in 8 byte chunks
        tmp_grad_data = []
        with zopen("131.0", mode="rb") as file:
            binary = file.read()
        tmp_grad_data.extend(
            struct.unpack("d", binary[ii * 8 : (ii + 1) * 8])[0]
            for ii in range(len(binary) // 8)
        )
        grad = [
            [
                float(tmp_grad_data[ii * 3]),
                float(tmp_grad_data[ii * 3 + 1]),
                float(tmp_grad_data[ii * 3 + 2]),
            ]
            for ii in range(len(tmp_grad_data) // 3)
        ]

        # Ensure that the scratch values match the correct values from the
        # output file but with higher precision
        if qc_output["pcm_gradients"] is not None:
            gradient = qc_output["pcm_gradients"][0]
        else:
            gradient = qc_output["gradients"][0]
        for ii, subgrad in enumerate(grad):
            for jj, val in enumerate(subgrad):
                if abs(gradient[ii, jj] - val) > 1e-6:
                    raise ValueError(
                        "Difference between gradient value in scratch file vs. output file should not be this large."
                    )
                gradient[ii, jj] = val

        # Convert gradient to force + deal with units
        results["forces"] = gradient * (-units.Hartree / units.Bohr)

    elif job_type == "freq":
        # Read Hessian scratch file in 8 byte chunks
        tmp_hess_data = []
        with zopen("132.0", mode="rb") as file:
            binary = file.read()
        tmp_hess_data.extend(
            struct.unpack("d", binary[ii * 8 : (ii + 1) * 8])[0]
            for ii in range(len(binary) // 8)
        )

        results["hessian"] = np.reshape(
            np.array(tmp_hess_data),
            (len(qc_output["species"]) * 3, len(qc_output["species"]) * 3),
        )
        results["enthalpy"] = qc_output["total_enthalpy"] * (units.kcal / units.mol)
        results["entropy"] = qc_output["total_entropy"] * (
            0.001 * units.kcal / units.mol
        )

    # Read orbital coefficients scratch file in 8 byte chunks
    prev_orbital_coeffs = []
    with zopen("53.0", mode="rb") as file:
        binary = file.read()
    prev_orbital_coeffs.extend(
        struct.unpack("d", binary[ii * 8 : (ii + 1) * 8])[0]
        for ii in range(len(binary) // 8)
    )

    return results, prev_orbital_coeffs

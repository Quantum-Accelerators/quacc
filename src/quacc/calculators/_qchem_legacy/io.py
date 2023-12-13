"""Q-Chem calculator IO."""
from __future__ import annotations

import struct
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np
from ase import units
from emmet.core.tasks import _parse_custodian
from monty.io import zopen
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.io.qchem.inputs import QCInput
from pymatgen.io.qchem.outputs import QCOutput
from pymatgen.io.qchem.sets import QChemDictSet

if TYPE_CHECKING:
    from typing import Any, Literal

    from ase.atoms import Atoms

    from quacc.calculators._qchem_legacy.qchem import Results


def write_qchem(
    atoms: Atoms,
    directory: Path | str = ".",
    charge: int = 0,
    spin_multiplicity: int = 1,
    basis_set: str = "def2-tzvpd",
    job_type: Literal["sp", "force", "opt", "freq"] = "force",
    scf_algorithm: str = "diis",
    qchem_input_params: dict[str, Any] | None = None,
    prev_orbital_coeffs: list[float] | None = None,
) -> None:
    """
    Write the Q-Chem input files.

    Parameters
    ----------
    atoms
        The Atoms object to be used for the calculation.
    directory
        The directory in which to write the Q-Chem input files.
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

    Returns
    -------
    None
    """
    directory = Path(directory)

    if prev_orbital_coeffs is not None:
        with Path(directory / "53.0").open(mode="wb") as file:
            for val in prev_orbital_coeffs:
                data = struct.pack("d", val)
                file.write(data)
        if "overwrite_inputs" not in qchem_input_params:
            qchem_input_params["overwrite_inputs"] = {}
        if "rem" not in qchem_input_params["overwrite_inputs"]:
            qchem_input_params["overwrite_inputs"]["rem"] = {}
        if "scf_guess" not in qchem_input_params["overwrite_inputs"]["rem"]:
            qchem_input_params["overwrite_inputs"]["rem"]["scf_guess"] = "read"

    atoms.charge = charge
    atoms.spin_multiplicity = spin_multiplicity
    mol = AseAtomsAdaptor.get_molecule(atoms)
    QChemDictSet(
        mol, job_type, basis_set, scf_algorithm, qchem_version=6, **qchem_input_params
    ).write(directory / "mol.qin")


def read_qchem(directory: Path | str = ".") -> tuple[Results, list[float]]:
    """
    Read Q-Chem log files.

    Parameters
    ----------
    directory
        The directory in which the Q-Chem calculation was run.

    Returns
    -------
    tuple[Results, list[float]]
        The results of the calculation and the orbital coefficients from a previous
        calculation.
    """
    directory = Path(directory)
    qc_input = QCInput.from_file(directory / "mol.qin").as_dict()
    qc_output = QCOutput(directory / "mol.qout").data

    results: Results = {
        "energy": qc_output["final_energy"] * units.Hartree,
        "qc_output": qc_output,
        "qc_input": qc_input,
        "custodian": _parse_custodian(directory),
    }

    # Read the gradient scratch file in 8 byte chunks
    grad_scratch = directory / "131.0"
    if grad_scratch.exists() and grad_scratch.stat().st_size > 0:
        tmp_grad_data = []
        with zopen(grad_scratch, mode="rb") as file:
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

    # Read Hessian scratch file in 8 byte chunks
    hessian_scratch = directory / "132.0"
    if hessian_scratch.exists() and hessian_scratch.stat().st_size > 0:
        tmp_hess_data = []
        with zopen(hessian_scratch, mode="rb") as file:
            binary = file.read()
        tmp_hess_data.extend(
            struct.unpack("d", binary[ii * 8 : (ii + 1) * 8])[0]
            for ii in range(len(binary) // 8)
        )
        reshaped_hess = np.reshape(
            np.array(tmp_hess_data),
            (len(qc_output["species"]) * 3, len(qc_output["species"]) * 3),
        )
        results["hessian"] = reshaped_hess * (units.Hartree / units.Bohr**2)
        results["enthalpy"] = qc_output["total_enthalpy"] * (units.kcal / units.mol)
        results["entropy"] = qc_output["total_entropy"] * (
            0.001 * units.kcal / units.mol
        )

    # Read orbital coefficients scratch file in 8 byte chunks
    prev_orbital_coeffs = []
    with zopen(directory / "53.0", mode="rb") as file:
        binary = file.read()
    prev_orbital_coeffs.extend(
        struct.unpack("d", binary[ii * 8 : (ii + 1) * 8])[0]
        for ii in range(len(binary) // 8)
    )

    return results, prev_orbital_coeffs

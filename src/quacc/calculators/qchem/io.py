from __future__ import annotations

import struct
from copy import deepcopy
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np
from ase import Atoms, units
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

    class Results(TypedDict):
        energy: float  # eV
        forces: NDArray  # Nx3, eV/Å
        hessian: NDArray  # Nx3x3, native Q-Chem units. TODO: convert to ASE units
        enthalpy: float  # eV
        entropy: float  # eV/K
        qc_output: dict[
            str, Any
        ]  # output from `pymatgen.io.qchem.outputs.QCOutput.data`
        qc_input: dict[
            str, Any
        ]  # input from `pymatgen.io.qchem.inputs.QCInput.as_dict()`
        custodian: dict[str, Any]  # custodian.json file metadata


def write_qchem(
    atoms: Atoms,
    directory: Path | str,
    input_filepath: Path | str,
    parameters: dict[str, Any],
) -> None:
    """
    Write the Q-Chem input files.

    Parameters
    ----------
    atoms
        The ASE atoms object.
    directory
        The directory to write the Q-Chem files to.
    input_filepath
        The path to the input file to write.
    paramters
        The Q-Chem parameters to use, formatted as a dictionary.

    Returns
    -------
    None
    """

    charge = parameters.get("charge", 0)
    spin_multiplicity = parameters.get("spin_multiplicity", 1)
    method = parameters.get("method")
    basis_set = parameters.get("basis_set", "def2-tzvpd")
    job_type = parameters.get("job_type", "force")
    scf_algorithm = parameters.get("scf_algorithm", "diis")
    qchem_input_params = parameters.get("qchem_input_params", {})

    atoms = deepcopy(atoms)
    atoms.charge = charge
    atoms.spin_multiplicity = spin_multiplicity
    mol = AseAtomsAdaptor.get_molecule(atoms)

    if self.prev_orbital_coeffs is not None:
        with Path(directory, "53.0").open(mode="wb") as file:
            for val in self.prev_orbital_coeffs:
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
    qcin.write(input_filepath)


def read_qchem(
    directory: Path | str,
    input_filepath: Path | str,
    ouptut_filepath: Path | str,
    job_type: Literal["opt", "sp", "freq", "force"],
) -> Results:
    """
    Read the Q-Chem output files.

    Parameters
    ----------
    directory
        The path to the directory to read the Q-Chem results from.
    input_filepath
        The path to the input file to read.
    output_filepath
        The path to the output file to read.
    job_type
        The type of Q-Chem job that was run.

    Returns
    -------
    Results
        The Q-Chem results, formatted as a dictionary.
    """
    qc_input = QCInput.from_file(str(input_filepath)).as_dict()
    qc_output = QCOutput(str(ouptut_filepath)).data

    results = {}
    results["qc_output"] = qc_output
    results["qc_input"] = qc_input
    results["custodian"] = _parse_custodian(directory)
    results["energy"] = qc_output["final_energy"] * units.Hartree

    if job_type in ["force", "opt"]:
        tmp_grad_data = []

        # Read the gradient scratch file in 8 byte chunks
        with zopen(directory / "131.0", mode="rb") as file:
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

    # Read orbital coefficients scratch file in 8 byte chunks
    self.prev_orbital_coeffs = []
    with zopen(directory / "53.0", mode="rb") as file:
        binary = file.read()
    self.prev_orbital_coeffs.extend(
        struct.unpack("d", binary[ii * 8 : (ii + 1) * 8])[0]
        for ii in range(len(binary) // 8)
    )

    if job_type == "freq":
        # Read Hessian scratch file in 8 byte chunks
        tmp_hess_data = []
        with zopen(directory / "132.0", mode="rb") as file:
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

"""Q-Chem calculator IO."""

from __future__ import annotations

import struct
from pathlib import Path
from typing import TYPE_CHECKING

from ase import units
from emmet.core.tasks import _parse_custodian
from pymatgen.io.qchem.inputs import QCInput
from pymatgen.io.qchem.outputs import (
    QCOutput,
    gradient_parser,
    hessian_parser,
    orbital_coeffs_parser,
)

if TYPE_CHECKING:
    from quacc.calculators.qchem.qchem import Results


def write_qchem(
    qc_input: QCInput,
    directory: Path | str = ".",
    prev_orbital_coeffs: list[float] | None = None,
) -> None:
    """
    Write the Q-Chem input files.

    Parameters
    ----------
    qc_input
        The QCInput object.
    directory
        The directory in which to write the files.
    prev_orbital_coeffs
        The orbital coefficients from a previous calculation.

    Returns
    -------
    None
    """
    directory = Path(directory)

    if prev_orbital_coeffs:
        with Path(directory / "53.0").open(mode="wb") as file:
            for val in prev_orbital_coeffs:
                data = struct.pack("d", val)
                file.write(data)

    qc_input.write_file(directory / "mol.qin")


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

    qc_input = QCInput.from_file(directory / "mol.qin")
    qc_output = QCOutput(directory / "mol.qout").data

    results: Results = {
        "energy": qc_output["final_energy"] * units.Hartree,
        "qc_output": qc_output,
        "qc_input": qc_input.as_dict(),
        "custodian": _parse_custodian(directory),
    }

    # Parse thermo properties
    if "total_enthalpy" in qc_output:
        results["enthalpy"] = qc_output["total_enthalpy"] * (units.kcal / units.mol)
    if "total_entropy" in qc_output:
        results["entropy"] = qc_output["total_entropy"] * (
            0.001 * units.kcal / units.mol
        )

    # Read the gradient scratch file in 8 byte chunks
    grad_scratch = directory / "131.0"
    if grad_scratch.exists() and grad_scratch.stat().st_size > 0:
        gradient = gradient_parser(directory / "131.0")

        results["forces"] = gradient * (-units.Hartree / units.Bohr)

    # Read Hessian scratch file in 8 byte chunks
    hessian_scratch = directory / "132.0"
    if hessian_scratch.exists() and hessian_scratch.stat().st_size > 0:
        reshaped_hess = hessian_parser(
            hessian_scratch, n_atoms=len(qc_output["species"])
        )
        results["hessian"] = reshaped_hess * (units.Hartree / units.Bohr**2)

    # Read orbital coefficients scratch file in 8 byte chunks
    orb_scratch = directory / "53.0"
    prev_orbital_coeffs = None
    if orb_scratch.exists() and orb_scratch.stat().st_size > 0:
        prev_orbital_coeffs = orbital_coeffs_parser(orb_scratch).tolist()

    return results, prev_orbital_coeffs

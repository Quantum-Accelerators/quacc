"""A Q-Chem calculator built on Pymatgen and Custodian functionality"""
from __future__ import annotations

import inspect
import logging
import struct
from copy import deepcopy
from pathlib import Path

import numpy as np
from ase import Atoms, units
from ase.calculators.calculator import FileIOCalculator
from emmet.core.tasks import _parse_custodian
from monty.io import zopen
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.io.qchem.inputs import QCInput
from pymatgen.io.qchem.outputs import QCOutput
from pymatgen.io.qchem.sets import QChemDictSet

from quacc.custodian import qchem as custodian_qchem

logger = logging.getLogger(__name__)


# TODO: Would be more consistent to convert the Hessian to ASE units if it's being
# stored in the first-class results.
class QChem(FileIOCalculator):
    """

    Parameters
    ----------
    atoms
        The Atoms object to be used for the calculation.
    cores
        Number of cores to use for the Q-Chem calculation.
    charge
        The total charge of the molecular system.
    spin_multiplicity
        The spin multiplicity of the molecular system.
    qchem_input_params
        Dictionary of Q-Chem input parameters to be passed to
        pymatgen.io.qchem.sets.DictSet.
    **fileiocalculator_kwargs
        Additional arguments to be passed to
        ase.calculators.calculator.FileIOCalculator.

    Returns
    -------
    Atoms
        The ASE Atoms object with attached Q-Chem calculator. In calc.results,
        the following properties are set:

        - energy: electronic energy in eV
        - forces: forces in eV/A
        - hessian: Hessian in native Q-Chem units
        - entropy: total enthalpy in eV
        - enthalpy: total entropy in eV/K
        - qc_output: Output from `pymatgen.io.qchem.outputs.QCOutput.data`
        - qc_input: Input from `pymatgen.io.qchem.inputs.QCInput.as_dict()`
        - custodian: custodian.json file metadata
    """

    implemented_properties = [
        "energy",
        "forces",
        "hessian",
        "enthalpy",
        "entropy",
        "qc_output",
        "qc_input",
        "custodian",
    ]  # noqa: RUF012

    def __init__(
        self,
        atoms: Atoms,
        charge: int = 0,
        spin_multiplicity: int = 1,
        method: str | None = None,
        basis_set: str = "def2-tzvpd",
        job_type: str = "force",
        scf_algorithm: str = "diis",
        cores: int = 1,
        qchem_input_params: dict | None = None,
        **fileiocalculator_kwargs,
    ):
        # Assign variables to self
        self.atoms = atoms
        self.cores = cores
        self.charge = charge
        self.spin_multiplicity = spin_multiplicity
        self.job_type = job_type
        self.basis_set = basis_set
        self.scf_algorithm = scf_algorithm
        self.qchem_input_params = qchem_input_params or {}
        self.fileiocalculator_kwargs = fileiocalculator_kwargs

        # Sanity checks
        if "directory" in self.fileiocalculator_kwargs:
            raise NotImplementedError("The directory kwarg is not supported.")

        if "overwrite_inputs" not in self.qchem_input_params:
            self.qchem_input_params["overwrite_inputs"] = {}

        if self.qchem_input_params.get("smd_solvent") and self.qchem_input_params.get(
            "pcm_dielectric"
        ):
            raise ValueError("PCM and SMD cannot be employed simultaneously.")

        if "rem" not in self.qchem_input_params["overwrite_inputs"]:
            self.qchem_input_params["overwrite_inputs"]["rem"] = {}
        if (
            method
            and "method" not in self.qchem_input_params["overwrite_inputs"]["rem"]
        ):
            self.qchem_input_params["overwrite_inputs"]["rem"]["method"] = method

        # We will save the parameters that have been passed to the Q-Chem
        # calculator via FileIOCalculator's self.default_parameters
        self.default_parameters = {
            "cores": self.cores,
            "charge": self.charge,
            "spin_multiplicity": self.spin_multiplicity,
            "scf_algorithm": self.scf_algorithm,
            "basis_set": self.basis_set,
        }

        if method:
            self.default_parameters["method"] = method

        # We also want to save the contents of self.qchem_input_params. However,
        # the overwrite_inputs key will have a corresponding value which is
        # either an empty dictionary or a nested dict of dicts, requiring a bit
        # of careful unwrapping.
        for key in self.qchem_input_params:
            if key == "overwrite_inputs":
                for subkey in self.qchem_input_params[key]:
                    for subsubkey in self.qchem_input_params[key][subkey]:
                        self.default_parameters[
                            f"overwrite_{subkey}_{subsubkey}"
                        ] = self.qchem_input_params[key][subkey][subsubkey]
            else:
                self.default_parameters[key] = self.qchem_input_params[key]

        # Get Q-Chem executable command
        self.command = self._manage_environment()

        # Instantiate previous orbital coefficients
        self.prev_orbital_coeffs = None

        # Instantiate the calculator
        FileIOCalculator.__init__(
            self,
            restart=None,
            ignore_bad_restart_file=FileIOCalculator._deprecated,
            label=None,
            atoms=self.atoms,
            **self.fileiocalculator_kwargs,
        )

    def _manage_environment(self) -> str:
        """
        Manage the environment for the Q-Chem calculator.

        Returns
        -------
        str
            The command flag to pass to the Q-Chem calculator.
        """

        # Return the command flag
        run_qchem_custodian_file = Path(inspect.getfile(custodian_qchem)).resolve()
        return f"python {run_qchem_custodian_file} {self.cores}"

    def write_input(self, atoms, properties=None, system_changes=None):
        FileIOCalculator.write_input(self, atoms, properties, system_changes)
        atoms = deepcopy(atoms)
        atoms.charge = self.charge
        atoms.spin_multiplicity = self.spin_multiplicity
        mol = AseAtomsAdaptor.get_molecule(atoms)
        if self.prev_orbital_coeffs is not None:
            with Path("53.0").open(mode="wb") as file:
                for val in self.prev_orbital_coeffs:
                    data = struct.pack("d", val)
                    file.write(data)
            if "overwrite_inputs" not in self.qchem_input_params:
                self.qchem_input_params["overwrite_inputs"] = {}
            if "rem" not in self.qchem_input_params["overwrite_inputs"]:
                self.qchem_input_params["overwrite_inputs"]["rem"] = {}
            if "scf_guess" not in self.qchem_input_params["overwrite_inputs"]["rem"]:
                self.qchem_input_params["overwrite_inputs"]["rem"]["scf_guess"] = "read"
        qcin = QChemDictSet(
            mol,
            self.job_type,
            self.basis_set,
            self.scf_algorithm,
            qchem_version=6,
            **self.qchem_input_params,
        )
        qcin.write("mol.qin")

    def read_results(self):
        qc_input = QCInput.from_file("mol.qin").as_dict()
        qc_output = QCOutput("mol.qout").data
        self.results["qc_output"] = qc_output
        self.results["qc_input"] = qc_input
        self.results["custodian"] = _parse_custodian(Path.cwd())

        self.results["energy"] = qc_output["final_energy"] * units.Hartree

        if self.job_type in ["force", "opt"]:
            tmp_grad_data = []
            # Read the gradient scratch file in 8 byte chunks
            with zopen("131.0", mode="rb") as file:
                binary = file.read()
            tmp_grad_data.extend(
                struct.unpack("d", binary[ii * 8 : (ii + 1) * 8])[0]
                for ii in range(int(len(binary) / 8))
            )
            grad = [
                [
                    float(tmp_grad_data[ii * 3]),
                    float(tmp_grad_data[ii * 3 + 1]),
                    float(tmp_grad_data[ii * 3 + 2]),
                ]
                for ii in range(int(len(tmp_grad_data) / 3))
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
            self.results["forces"] = gradient * (-units.Hartree / units.Bohr)
        else:
            self.results["forces"] = None

        # Read orbital coefficients scratch file in 8 byte chunks
        self.prev_orbital_coeffs = []
        with zopen("53.0", mode="rb") as file:
            binary = file.read()
        self.prev_orbital_coeffs.extend(
            struct.unpack("d", binary[ii * 8 : (ii + 1) * 8])[0]
            for ii in range(int(len(binary) / 8))
        )

        if self.job_type == "freq":
            tmp_hess_data = []
            # Read Hessian scratch file in 8 byte chunks
            with zopen("132.0", mode="rb") as file:
                binary = file.read()
            tmp_hess_data.extend(
                struct.unpack("d", binary[ii * 8 : (ii + 1) * 8])[0]
                for ii in range(int(len(binary) / 8))
            )
            self.results["hessian"] = np.reshape(
                np.array(tmp_hess_data),
                (len(qc_output["species"]) * 3, len(qc_output["species"]) * 3),
            )
            self.results["enthalpy"] = qc_output["total_enthalpy"] * (
                units.kcal / units.mol
            )
            self.results["entropy"] = qc_output["total_entropy"] * (
                0.001 * units.kcal / units.mol
            )
        else:
            self.results["hessian"] = None

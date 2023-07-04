"""
A Q-Chem calculator built on Pymatgen and Custodian functionality
"""
from __future__ import annotations

import inspect
import os
import struct

from ase import Atoms, units
from ase.calculators.calculator import FileIOCalculator
from monty.io import zopen
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.io.qchem.outputs import QCOutput
from pymatgen.io.qchem.sets import ForceSet

from quacc.custodian import qchem as custodian_qchem
from quacc.util.atoms import check_charge_and_spin


class QChem(FileIOCalculator):
    """

    Parameters
    ----------
    input_atoms
        The Atoms object to be used for the calculation.
    cores
        Number of cores to use for the Q-Chem calculation.
    charge
        The total charge of the molecular system.
        Effectively defaults to zero.
    spin_multiplicity
        The spin multiplicity of the molecular system.
        Effectively defaults to the lowest spin state given the molecular structure and charge.
    qchem_input_params
        Dictionary of Q-Chem input parameters to be passed to pymatgen.io.qchem.sets.ForceSet.
    **fileiocalculator_kwargs
        Additional arguments to be passed to ase.calculators.calculator.FileIOCalculator.

    Returns
    -------
    Atoms
        The ASE Atoms object with attached Q-Chem calculator.
    """

    implemented_properties = ["energy", "forces"]

    def __init__(
        self,
        input_atoms: Atoms,
        charge: None | int = None,
        spin_multiplicity: None | int = None,
        method: str | None = None,
        cores: None | int = None,
        qchem_input_params: dict | None = None,
        **fileiocalculator_kwargs,
    ):
        # Assign variables to self
        self.input_atoms = input_atoms
        self.cores = cores
        self.charge = charge
        self.spin_multiplicity = spin_multiplicity
        self.qchem_input_params = qchem_input_params or {}
        self.fileiocalculator_kwargs = fileiocalculator_kwargs

        # Sanity checks
        if "overwrite_inputs" not in self.qchem_input_params:
            self.qchem_input_params["overwrite_inputs"] = {}

        if self.charge is None and self.spin_multiplicity is not None:
            raise ValueError(
                "If setting spin_multiplicity, must also specify charge! Exiting..."
            )

        if self.qchem_input_params.get("smd_solvent") and self.qchem_input_params.get(
            "pcm_dielectric"
        ):
            raise ValueError(
                "PCM and SMD cannot be employed simultaneously! Exiting..."
            )

        if "rem" not in self.qchem_input_params["overwrite_inputs"]:
            self.qchem_input_paramsoverwrite_inputs["rem"] = {}
        if (
            method
            and "method" not in self.qchem_input_params["overwrite_inputs"]["rem"]
        ):
            self.qchem_input_params["overwrite_inputs"]["rem"]["method"] = method

        # We will save the parameters that have been passed to the Q-Chem calculator via FileIOCalculator's
        # self.default_parameters
        self.default_parameters = {
            "cores": self.cores,
            "charge": self.charge,
            "spin_multiplicity": self.spin_multiplicity,
        }

        # We also want to save the contents of self.qchem_input_params. However, the overwrite_inputs
        # key will have a corresponding value which is either an empty dictionary or a nested dict of
        # dicts, requiring a bit of careful unwrapping.
        for key in self.qchem_input_params:
            if key == "overwrite_inputs":
                for subkey in self.qchem_input_params[key]:
                    for subsubkey in self.qchem_input_params[key][subkey]:
                        self.default_parameters[
                            f"overwrite_{subkey}_{subsubkey}"
                        ] = self.qchem_input_params[key][subkey][subsubkey]
            else:
                self.default_parameters[key] = self.qchem_input_params[key]

        charge, spin_multiplicity = check_charge_and_spin(
            input_atoms, self.charge, self.spin_multiplicity
        )
        self.charge = charge
        self.spin_multiplicity = spin_multiplicity

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
            atoms=self.input_atoms,
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
        run_qchem_custodian_file = os.path.abspath(inspect.getfile(custodian_qchem))
        if self.cores is not None:
            return f"python {run_qchem_custodian_file} {self.cores}"
        else:
            return f"python {run_qchem_custodian_file}"

    def write_input(self, atoms, properties=None, system_changes=None):
        FileIOCalculator.write_input(self, atoms, properties, system_changes)
        atoms.charge = self.charge
        atoms.spin_multiplicity = self.spin_multiplicity
        mol = AseAtomsAdaptor.get_molecule(atoms)
        if self.prev_orbital_coeffs is not None:
            with open("53.0", mode="wb") as file:
                for val in self.prev_orbital_coeffs:
                    data = struct.pack("d", val)
                    file.write(data)
            if "overwrite_inputs" not in self.qchem_input_params:
                self.qchem_input_params["overwrite_inputs"] = {}
            if "rem" not in self.qchem_input_params["overwrite_inputs"]:
                self.qchem_input_params["overwrite_inputs"]["rem"] = {}
            if "scf_guess" not in self.qchem_input_params["overwrite_inputs"]["rem"]:
                self.qchem_input_params["overwrite_inputs"]["rem"]["scf_guess"] = "read"
        qcin = ForceSet(mol, qchem_version=6, **self.qchem_input_params)
        qcin.write("mol.qin")

    def read_results(self):
        data = QCOutput("mol.qout").data
        self.results["energy"] = data["final_energy"] * units.Hartree
        tmp_grad_data = []
        # Read the gradient scratch file in 8 byte chunks
        with zopen("131.0", mode="rb") as file:
            binary = file.read()
            for ii in range(int(len(binary) / 8)):
                tmp_grad_data.append(
                    struct.unpack("d", binary[ii * 8 : (ii + 1) * 8])[0]
                )
        # Reshape the gradient into N x 3
        grad = []
        for ii in range(int(len(tmp_grad_data) / 3)):
            grad.append(
                [
                    float(tmp_grad_data[ii * 3]),
                    float(tmp_grad_data[ii * 3 + 1]),
                    float(tmp_grad_data[ii * 3 + 2]),
                ]
            )
        # Ensure that the scratch values match the correct values from the output file
        # but with higher precision
        if data["pcm_gradients"] is not None:
            gradient = data["pcm_gradients"][0]
        else:
            gradient = data["gradients"][0]
        for ii, subgrad in enumerate(grad):
            for jj, val in enumerate(subgrad):
                if abs(gradient[ii, jj] - val) > 1e-6:
                    raise ValueError(
                        "Difference between gradient value in scratch file vs. output file should not be this large! Exiting..."
                    )
                gradient[ii, jj] = val
        # Convert gradient to force + deal with units
        self.results["forces"] = gradient * (-units.Hartree / units.Bohr)
        self.prev_orbital_coeffs = []
        # Read orbital coefficients scratch file in 8 byte chunks
        with zopen("53.0", mode="rb") as file:
            binary = file.read()
            for ii in range(int(len(binary) / 8)):
                self.prev_orbital_coeffs.append(
                    struct.unpack("d", binary[ii * 8 : (ii + 1) * 8])[0]
                )

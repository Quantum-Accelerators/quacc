"""
A wrapper around ASE's Q-Chem calculator that makes it better suited for high-throughput DFT.
"""
from __future__ import annotations

import inspect
import os
import struct
import warnings

import ase.units
import numpy as np
from ase import Atoms
from ase.calculators.calculator import FileIOCalculator
from monty.io import zopen
from pymatgen.core import Molecule
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.io.qchem.inputs import QCInput
from pymatgen.io.qchem.outputs import QCOutput
from pymatgen.io.qchem.sets import ForceSet

from quacc import SETTINGS
from quacc.custodian import qchem as custodian_qchem


class QChem(FileIOCalculator):
    """

    Parameters
    ----------
    input_atoms
        The Atoms object to be used for the calculation.
    use_custodian
        Whether to use Custodian to run Q-Chem.
        Default is True in settings.
    **kwargs
        Additional arguments to be passed to the Q-Chem calculator. Takes all valid
        ASE calculator arguments, in addition to those custom to Quacc.

    Returns
    -------
    Atoms
        The ASE Atoms object with attached Q-Chem calculator.
    """

    def __init__(
        self,
        input_atoms: Atoms,
        charge: None | int = None,
        spin_multiplicity: None | int = None,
        qchem_input_params: dict = {},
        use_custodian: bool = SETTINGS.QChem_CUSTODIAN,
        **kwargs,
    ):
        # Assign variables to self
        self.input_atoms = input_atoms
        self.charge = charge
        self.spin_multiplicity = spin_multiplicity
        self.qchem_input_params = qchem_input_params
        self.use_custodian = use_custodian
        self.kwargs = kwargs

        # Instantiate the calculator
        FileIOCalculator.__init__(
            self,
            restart=None,
            ignore_bad_restart_file=FileIOCalculator._deprecated,
            label=None,
            atoms=self.input_atoms,
            **self.kwargs,
        )

        # Get Q-Chem executable command
        command = self._manage_environment()

    def _manage_environment(self) -> str:
        """
        Manage the environment for the Q-Chem calculator.

        Returns
        -------
        str
            The command flag to pass to the Q-Chem calculator.
        """

        # Check if Custodian should be used
        if self.use_custodian:
            # Return the command flag
            run_qchem_custodian_file = os.path.abspath(inspect.getfile(custodian_qchem))
            return f"python {run_qchem_custodian_file}"
        else:
            warnings.warn("Can only run Q-Chem via Custodian currently!")
            return None

    def write_input(self, atoms, properties=None, system_changes=None):
        FileIOCalculator.write_input(self, atoms, properties, system_changes)
        mol = AseAtomsAdaptor.get_molecule(atoms)
        if self.charge is not None:
            if self.spin_multiplicity is not None:
                mol.set_charge_and_spin(
                    charge=self.charge, spin_multiplicity=self.spin_multiplicity
                )
            else:
                mol.set_charge_and_spin(charge=self.charge)
        qcin = ForceSet(mol, **self.qchem_input_params)
        qcin.write("mol.qin")

    def read_results(self):
        data = QCOutput("mol.qout").data
        self.results["energy"] = data["final_energy"] * ase.units.Hartree
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
        gradient = data["gradients"][0]
        for ii, subgrad in enumerate(grad):
            for jj, val in enumerate(subgrad):
                assert abs(gradient[ii, jj] - val) < 1e-7
                gradient[ii, jj] = val
        self.results["forces"] = gradient * (-ase.units.Hartree / ase.units.Bohr)

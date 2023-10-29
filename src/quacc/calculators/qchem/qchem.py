"""A Q-Chem calculator built on Pymatgen and Custodian functionality"""
from __future__ import annotations

import inspect
from pathlib import Path
from typing import TYPE_CHECKING

from ase.calculators.calculator import FileIOCalculator

from quacc.calculators.qchem import custodian
from quacc.calculators.qchem.io import read_qchem, write_qchem
from quacc.calculators.qchem.params import cleanup_params, get_default_params

if TYPE_CHECKING:
    from typing import Any, ClassVar, Literal

    from ase import Atoms

    from quacc.calculators.qchem.io import Results


class QChem(FileIOCalculator):
    implemented_properties: ClassVar[list[str]] = [
        "energy",
        "forces",
        "hessian",
        "enthalpy",
        "entropy",
        "qc_output",
        "qc_input",
        "custodian",
    ]
    results: ClassVar[Results] = {}

    def __init__(
        self,
        atoms: Atoms,
        charge: int = 0,
        spin_multiplicity: int = 1,
        method: str | None = None,
        basis_set: str = "def2-tzvpd",
        job_type: Literal["sp", "force", "opt", "freq"] = "force",
        scf_algorithm: str = "diis",
        cores: int = 1,
        qchem_input_params: dict[str, Any] | None = None,
        **fileiocalculator_kwargs,
    ) -> None:
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
            `pymatgen.io.qchem.sets.DictSet`.
        **fileiocalculator_kwargs
            Additional arguments to be passed to
            `ase.calculators.calculator.FileIOCalculator`.

        Returns
        -------
        Atoms
            The ASE Atoms object with attached Q-Chem calculator.
        """

        # Assign variables to self
        self.atoms = atoms
        self.cores = cores
        self.charge = charge
        self.spin_multiplicity = spin_multiplicity
        self.job_type = job_type
        self.basis_set = basis_set
        self.scf_algorithm = scf_algorithm
        self.qchem_input_params = (
            qchem_input_params if qchem_input_params is not None else {}
        )
        self.fileiocalculator_kwargs = fileiocalculator_kwargs

        # Instantiate previous orbital coefficients
        self.prev_orbital_coeffs = None

        if "directory" in self.fileiocalculator_kwargs:
            raise NotImplementedError("The directory kwarg is not supported.")

        # Set parameters
        self.qchem_input_params = cleanup_params(self.qchem_input_params, method)
        self.default_parameters = get_default_params(
            cores,
            charge,
            spin_multiplicity,
            method,
            basis_set,
            scf_algorithm,
            qchem_input_params,
        )

        # Get Q-Chem executable command
        self.command = self._manage_environment()

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
            The command flag to run Q-Chem with Custodian.
        """

        qchem_custodian_script = Path(inspect.getfile(custodian)).resolve()
        return f"python {qchem_custodian_script} {self.cores}"

    def write_input(
        self,
        atoms: Atoms,
        properties: list[str] | None = None,
        system_changes: list[str] | None = None,
    ) -> None:
        """
        Write the Q-Chem input files.

        Parameters
        ----------
        atoms
            The Atoms object to be used for the calculation.
        properties
            List of properties to calculate.
        system_changes
            List of changes to the system since last calculation.

        Returns
        -------
        None
        """
        write_qchem(
            atoms,
            charge=self.charge,
            spin_multiplicity=self.spin_multiplicity,
            basis_set=self.basis_set,
            job_type=self.job_type,
            scf_algorithm=self.scf_algorithm,
            qchem_input_params=self.qchem_input_params,
            prev_orbital_coeffs=self.prev_orbital_coeffs,
            properties=properties,
            system_changes=system_changes,
        )

    def read_results(self) -> None:
        """
        Read the Q-Chem output files.

        Parameters
        ----------
        None

        Returns
        -------
        None
        """
        results, prev_orbital_coeffs = read_qchem(self.job_type)
        self.results = results
        self.prev_orbital_coeffs = prev_orbital_coeffs

"""A Q-Chem calculator built on Pymatgen and Custodian functionality."""
from __future__ import annotations

import inspect
from pathlib import Path
from typing import TYPE_CHECKING

from ase.calculators.calculator import FileIOCalculator

from quacc.calculators._qchem_legacy import qchem_custodian
from quacc.calculators._qchem_legacy.io import read_qchem, write_qchem

if TYPE_CHECKING:
    from typing import Any, ClassVar, Literal, TypedDict

    from ase.atoms import Atoms
    from numpy.typing import NDArray

    class Results(TypedDict, total=False):
        energy: float  # electronic energy in eV
        forces: NDArray  # forces in eV/A
        hessian: NDArray  # Hessian in eV/A^2/amu
        enthalpy: float  # total enthalpy in eV
        entropy: float  # total entropy in eV/K
        qc_output: dict[
            str, Any
        ]  # Output from `pymatgen.io.qchem.outputs.QCOutput.data`
        qc_input: dict[
            str, Any
        ]  # Input from `pymatgen.io.qchem.inputs.QCInput.as_dict()`
        custodian: dict[str, Any]  # custodian.json file metadata


class QChem(FileIOCalculator):
    """Custom Q-Chem calculator built on Pymatgen and Custodian."""

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
        Initialize the Q-Chem calculator.

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
        method
            The level of theory to use.
        basis_set
            The basis set.
        job_type
            The job type for the calculation.
        scf_algorithm
            The SCF algorithm to use
        cores
            The number of CPU cores to run on.
        qchem_input_params
            Dictionary of Q-Chem input parameters to be passed to
            `pymatgen.io.qchem.sets.DictSet`.
        **fileiocalculator_kwargs
            Additional arguments to be passed to
            `ase.calculators.calculator.FileIOCalculator`.

        Returns
        -------
        None
        """

        # Assign variables to self
        self.atoms = atoms
        self.charge = charge
        self.spin_multiplicity = spin_multiplicity
        self.method = method
        self.basis_set = basis_set
        self.job_type = job_type
        self.scf_algorithm = scf_algorithm
        self.cores = cores
        self.qchem_input_params = qchem_input_params or {}
        self.fileiocalculator_kwargs = fileiocalculator_kwargs

        # Instantiate previous orbital coefficients
        self._prev_orbital_coeffs = None

        if "directory" in self.fileiocalculator_kwargs:
            raise NotImplementedError("The directory kwarg is not supported.")

        # Clean up parameters
        self._cleanup_qchem_input_params()
        self._set_default_params()

        # Get Q-Chem executable command
        command = self._manage_environment()

        # Instantiate the calculator
        super().__init__(
            restart=None,
            ignore_bad_restart_file=FileIOCalculator._deprecated,
            label=None,
            atoms=self.atoms,
            command=command,
            profile=None,
            **self.fileiocalculator_kwargs,
        )

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
        FileIOCalculator.write_input(self, atoms, properties, system_changes)
        write_qchem(
            atoms,
            charge=self.charge,
            spin_multiplicity=self.spin_multiplicity,
            basis_set=self.basis_set,
            job_type=self.job_type,
            scf_algorithm=self.scf_algorithm,
            qchem_input_params=self.qchem_input_params,
            prev_orbital_coeffs=self._prev_orbital_coeffs,
        )

    def read_results(self) -> None:
        """
        Read the Q-Chem output files. Update the .results and ._prev_orbital_coeffs
        attributes.

        Returns
        -------
        None
        """
        results, _prev_orbital_coeffs = read_qchem()
        self.results = results
        self._prev_orbital_coeffs = _prev_orbital_coeffs

    def _manage_environment(self) -> str:
        """
        Return the command to run the Q-Chem calculator via Custodian.

        Returns
        -------
        str
            The command flag to run Q-Chem with Custodian.
        """

        qchem_custodian_script = Path(inspect.getfile(qchem_custodian)).resolve()
        return f"python {qchem_custodian_script} {self.cores}"

    def _cleanup_qchem_input_params(self) -> None:
        """
        Clean up q-chem input parameters for the Q-Chem calculator. Modifies
        self.qchem_input_params in place.

        Parameters
        ----------
        None

        Returns
        -------
        None
        """
        if "overwrite_inputs" not in self.qchem_input_params:
            self.qchem_input_params["overwrite_inputs"] = {}

        if self.qchem_input_params.get("smd_solvent") and self.qchem_input_params.get(
            "pcm_dielectric"
        ):
            raise ValueError("PCM and SMD cannot be employed simultaneously.")

        if "rem" not in self.qchem_input_params["overwrite_inputs"]:
            self.qchem_input_params["overwrite_inputs"]["rem"] = {}
        if (
            self.method
            and "method" not in self.qchem_input_params["overwrite_inputs"]["rem"]
        ):
            self.qchem_input_params["overwrite_inputs"]["rem"]["method"] = self.method

    def _set_default_params(self) -> None:
        """
        Store the parameters that have been passed to the Q-Chem calculator in
        FileIOCalculator's self.default_parameters.

        Parameters
        ----------
        None

        Returns
        -------
        None
        """
        self.default_parameters = {
            "cores": self.cores,
            "charge": self.charge,
            "spin_multiplicity": self.spin_multiplicity,
            "scf_algorithm": self.scf_algorithm,
            "basis_set": self.basis_set,
        }

        if self.method:
            self.default_parameters["method"] = self.method

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

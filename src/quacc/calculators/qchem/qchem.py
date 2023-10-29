"""
ASE calculator for Q-Chem
"""
from __future__ import annotations

import inspect
import multiprocessing
from pathlib import Path
from subprocess import check_call
from typing import TYPE_CHECKING

from ase.calculators.genericfileio import CalculatorTemplate, GenericFileIOCalculator

from quacc.calculators.qchem import custodian
from quacc.calculators.qchem.io import read_qchem, write_qchem

if TYPE_CHECKING:
    from typing import Any

    from ase.atoms import Atoms

    from quacc.calculators.qchem.io import Results

class QChemProfile:
    """
    Q-Chem profile
    """

    def __init__(self, cores: int | None = None) -> None:
        """
        Initialize the Q-Chem profile.

        Parameters
        ----------
        cores
            The number of cores to use for the Q-Chem calculation.
            Defaults to the number of cores on the machine.

        Returns
        -------
        None
        """
        self.cores = cores or multiprocessing.cpu_count()

    def run(
        self,
        directory: Path | str,
        output_filename: str,
    ) -> None:
        """
        Run the xTB calculation.

        Parameters
        ----------
        directory
            The directory where the calculation will be run.
        output_filename
            The name of the log file to write to in the directory.

        Returns
        -------
        None
        """
        run_qchem_custodian_file = Path(inspect.getfile(custodian)).resolve()
        cmd = ["python", str(run_qchem_custodian_file), str(self.cores)]

        with open(output_filename, "w") as fd:
            check_call(cmd, stdout=fd, cwd=directory)


class _QChemTemplate(CalculatorTemplate):
    """
    Q-Chem template
    """

    def __init__(self) -> None:
        """
        Initialize the Q-Chem template.

        Parameters
        ----------
        None

        Returns
        -------
        None
        """
        label = "mol"
        super().__init__(
            name=label,
            implemented_properties=[
                "energy",
                "forces",
                "hessian",
                "enthalpy",
                "entropy",
                "qc_output",
                "qc_input",
                "custodian",
            ],
        )

        self.input_file = f"{label}.qin"
        self.output_file = f"{label}.qout"

    def execute(self, directory: Path | str, profile: QChemProfile) -> None:
        """
        Run the Q-Chem executable.

        Parameters
        ----------
        directory
            The path to the directory to run the Q-Chem executable in.
        profile
            The Q-Chem profile to use.

        Returns
        -------
        None
        """
        profile.run(directory, self.output_file)

    def write_input(
        self,
        directory: Path | str,
        atoms: Atoms,
        parameters: dict[str, Any],
        properties: Any,  # skipcq: PYL-W0613
    ) -> None:
        """
        Write the Q-Chem input files.

        Parameters
        ----------
        directory
            The path to the directory to write the Q-Chem input files in.
        atoms
            The ASE atoms object to write.
        parameters
            The Q-Chem parameters to use, formatted as a dictionary.
        properties
            This is needed the base class and should not be explicitly specified.

        Returns
        -------
        None
        """

        write_qchem(
            atoms,
            directory,
            directory / self.input_file,
            parameters,
        )

    def read_results(self, directory: Path) -> Results:
        """
        Use cclib to read the results from the Q-Chem calculation.

        Parameters
        ----------
        directory
            The path to the directory to read the Q-Chem results from.

        Returns
        -------
        Results
            The Q-Chem results, formatted as a dictionary.
        """

        return read_qchem(
            directory, directory / self.input_file, directory / self.output_file
        )


class QChem(GenericFileIOCalculator):
    """
    Q-Chem calculator
    """

    def __init__(
        self,
        *,
        profile: QChemProfile | None = None,
        directory: Path | str = ".",
        **kwargs,
    ) -> None:
        """
        Initialize the Q-Chem calculator.

        Parameters
        ----------
        profile
            The Q-Chem profile to use.
        directory
            The path to the directory to run the Q-Chem executable in.
        **kwargs
            The Q-Chem parameters to use:

            atoms: Atoms
                The Atoms object to be used for the calculation.
            charge: int
                The total charge of the molecular system.
                Default: 0.
            spin_multiplicity: int
                The spin multiplicity of the molecular system.
                Default: 1.
            method: str
                The Q-Chem method to use.
            basis_set: str
                The Q-Chem basis set to use.
            job_type: Literal["opt", "sp", "freq", "force"])
                The Q-Chem job type to use.
            scf_algorithm: str
                The Q-Chem SCF algorithm to use.
            qchem_input_params: dict[str, Any]
                Dictionary of Q-Chem input parameters to be passed to
                `pymatgen.io.qchem.sets.DictSet`.

        Returns
        -------
        None
        """

        profile = profile or QChemProfile()
        directory = Path(directory).expanduser().resolve()

        super().__init__(
            template=_QChemTemplate(),
            profile=profile,
            directory=directory,
            parameters=kwargs,
        )

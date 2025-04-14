from __future__ import annotations

import re
from typing import TYPE_CHECKING

from ase.calculators.genericfileio import (
    BaseProfile,
    CalculatorTemplate,
    GenericFileIOCalculator,
    read_stdout,
)

from quacc.calculators.mrcc.io import read_mrcc_outputs, write_mrcc

if TYPE_CHECKING:
    from pathlib import Path
    from typing import Any

    from ase.atoms import Atoms
    from ase.config import Config

    from quacc.calculators.mrcc.io import MRCCEnergyInfo


def _get_version_from_mrcc_header(mrcc_header: str) -> str:
    """
    Get the version of MRCC from the header of the output file.

    Parameters
    ----------
    mrcc_header
        The header of the MRCC output file.

    Returns
    -------
    str
        The version of the MRCC executable normally in Month DD, YYYY.
    """

    match = re.search(r"Release date: (.*)$", mrcc_header, re.M)
    return match.group(1)


class MrccProfile(BaseProfile):
    """The MRCC calculator profile"""

    def version(self) -> str:
        """
        Gives the version of the MRCC executable.

        Returns
        -------
        str
            The version of the MRCC executable normally in Month DD, YYYY.
        """
        stdout = read_stdout([self.command, "does_not_exist"])
        return _get_version_from_mrcc_header(stdout)

    def get_calculator_command(self, inputfile: str) -> list[str]:
        """
        Get the command to run the MRCC calculation.

        Parameters
        ----------
        inputfile
            The input file to run the calculation.

        Returns
        -------
        list[str]
            The command to run the MRCC calculation.
        """
        return [inputfile]


class MrccTemplate(CalculatorTemplate):
    """
    The MRCC calculator template class.
    """

    _label = "mrcc"

    def __init__(self) -> None:
        """
        Initialize the MRCC calculator template class.

        Returns
        -------
        None
        """
        super().__init__(
            "mrcc",
            implemented_properties=[
                "energy",
                "scf_energy",
                "mp2_corr_energy",
                "ccsd_corr_energy",
                "ccsdt_corr_energy",
            ],
        )
        self.inputname = "MINP"
        self.outputname = f"{self._label}.out"
        self.errorname = f"{self._label}.err"

    def execute(self, directory: Path | str, profile: MrccProfile) -> None:
        """
        Execute the MRCC calculation.

        Parameters
        ----------
        directory
            The directory in which to run the calculation.
        profile
            The MRCCProfile class.

        Returns
        -------
        None
        """
        profile.run(
            directory, self.inputname, self.outputname, errorfile=self.errorname
        )

    def write_input(
        self,
        profile: MrccProfile,  # noqa: ARG002
        directory: Path | str,
        atoms: Atoms,
        parameters: dict[str, str],
        properties: dict[str, Any],  # noqa: ARG002
    ) -> None:
        """
        Write the MRCC input file.

        Parameters
        ----------
        profile
            The MRCCProfile class.
        directory
            The directory in which to write the input file.
        atoms
            The Atoms object.
        parameters
            The parameters for the calculation.
        properties
            The properties to calculate.

        Returns
        -------
        None
        """
        parameters = dict(parameters)

        kw = {"charge": 0, "mult": 1, "calc": "PBE", "basis": "def2-SVP"}
        kw |= parameters

        write_mrcc(directory / self.inputname, atoms, kw)

    def read_results(self, directory: Path | str) -> MRCCEnergyInfo:
        """
        Reads the MRCC output files.

        Parameters
        ----------
        directory
            The directory in which the calculation was run.

        Returns
        -------
        MRCCEnergyInfo
            Dictionary with the energy components. The keys are the following:
            - energy : float <-- Total energy which will not be computed in this function.
            - scf_energy : float <-- SCF energy.
            - mp2_corr_energy : float <-- MP2 correlation energy.
            - ccsd_corr_energy : float <-- CCSD correlation energy.
            - ccsdt_corr_energy : float <-- CCSD(T) correlation energy.
        """
        return read_mrcc_outputs(output_file_path=directory / self.outputname)

    def load_profile(self, cfg: Config, **kwargs) -> MrccProfile:
        return MrccProfile.from_config(cfg, self.name, **kwargs)


class MRCC(GenericFileIOCalculator):
    """
    Class for performing MRCC calculations.
    """

    def __init__(
        self, *, profile: MrccProfile = None, directory: str | Path = ".", **kwargs
    ) -> None:
        """
        Construct MRCC-calculator object.

        Parameters
        ----------
        profile: MrccProfile
            The MRCC profile to use.
        directory: str
            The directory in which to run the calculation.
        **kwargs
            The parameters for the MRCC calculation.

        Examples
        --------
        Use default values:

        >>> from quacc.calculators.mrcc.mrcc import MRCC, MrccProfile
        >>> from ase.build import molecule
        >>> from quacc import get_settings

        >>> calc = MRCC(
        ...     profile=MrccProfile(command=get_settings().MRCC_CMD),
        ...     charge=0,
        ...     mult=1,
        ...     basis="def2-SVP",
        ...     calc="PBE",
        ... )
        >>> h = molecule("H2")
        >>> h.set_calculator(calc)
        >>> h.get_total_energy()

        Returns
        -------
        None
        """

        super().__init__(
            template=MrccTemplate(),
            profile=profile,
            directory=directory,
            parameters=kwargs,
        )

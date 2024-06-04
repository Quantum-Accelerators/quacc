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

    from ase import Atoms

    from quacc.calculators.mrcc.io import ParamsInfo


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
        List[str]
            The command to run the MRCC calculation.
        """
        return [inputfile]


class MrccTemplate(CalculatorTemplate):
    _label = "mrcc"

    def __init__(self):
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
        self, *, directory: Path | str, atoms: Atoms, parameters: ParamsInfo, **kwargs
    ) -> None:
        """
        Write the MRCC input file.

        Parameters
        ----------
        profile
            The MRCC profile.
        directory
            The directory in which to write the input file.
        atoms
            The Atoms object.
        parameters
            The parameters for the calculation.

        Returns
        -------
        None
        """
        parameters = dict(parameters)

        mrccinput = {"calc": "PBE", "basis": "def2-SVP"}

        all_parameters = {
            "charge": 0,
            "mult": 1,
            "mrccinput": mrccinput,
            "mrccblocks": "",
        }
        all_parameters.update(parameters)

        write_mrcc(
            file_path=directory / self.inputname, atoms=atoms, parameters=all_parameters
        )

    def read_results(self, directory: Path | str):
        """
        Reads the MRCC output files.

        Parameters
        ----------
        directory
            The directory in which the calculation was run.

        Returns
        -------
        EnergyInfo
            Dictionary with the energy components. The keys are the following:
            - energy : float <-- Total energy which will not be computed in this function.
            - scf_energy : float <-- SCF energy.
            - mp2_corr_energy : float <-- MP2 correlation energy.
            - ccsd_corr_energy : float <-- CCSD correlation energy.
            - ccsdt_corr_energy : float <-- CCSD(T) correlation energy.

        """
        return read_mrcc_outputs(output_file_path=directory / self.outputname)

    def load_profile(self, cfg, **kwargs):
        return MrccProfile.from_config(cfg, self.name, **kwargs)


class MRCC(GenericFileIOCalculator):
    """
    Class for performing MRCC calculations.

    Example:

      calc = MRCC(charge=0, mult=1, mrccinput={'calc':PBE, 'basis':'def2-SVP'},
        mrccblocks=' ')
    """

    def __init__(self, *, profile=None, directory=".", **kwargs):
        """Construct MRCC-calculator object.

        Parameters
        ----------
        charge: int
            The charge of the system.
        mult: int
            The multiplicity of the system.
        mrccinput : dict[str,str]
            The input for the MRCC calculation. The keys are the MRCC settings with the values being the corresponding value for each setting.
        mrccblocks: str
            The MRCC blocks to be written that goes after mrccinput.
        profile: MrccProfile
            The MRCC profile to use.
        directory: str
            The directory in which to run the calculation.

        """

        super().__init__(
            template=MrccTemplate(),
            profile=profile,
            directory=directory,
            parameters=kwargs,
        )

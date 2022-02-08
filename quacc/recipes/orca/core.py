from dataclasses import dataclass
from typing import Any, Dict

from ase.atoms import Atoms
from ase.calculators.orca import ORCA
from jobflow import Maker, job

from quacc.schemas.cclib import summarize_run
from quacc.util.calc import run_calc


@dataclass
class StaticMaker(Maker):
    """
    Class to carry out a single-point calculation.

    Parameters
    ----------
    name
        Name of the job.
    orcasimpleinput
        ORCA simple input string.
    orcablock
        ORCA block input string.
    """

    name: str = "ORCA-Static"
    orcasimpleinput: str = "SP SlowConv NormalPrint"
    orcablock: str = ""

    @job
    def make(self, atoms: Atoms) -> Dict[str, Any]:
        """
        Make the run.

        Parameters
        ----------
        atoms
            .Atoms object

        Returns
        -------
        Dict
            Summary of the run.
        """
        atoms.calc = ORCA(
            orcasimpleinput=self.orcasimpleinput,
            orcablock=self.orcablock,
        )
        atoms = run_calc(atoms)
        summary = summarize_run(atoms, ".out", additional_fields={"name": self.name})

        return summary


@dataclass
class RelaxMaker(Maker):
    """
    Class to carry out a single-point calculation.

    Parameters
    ----------
    name
        Name of the job.
    orcasimpleinput
        ORCA simple input string.
    orcablock
        ORCA block input string.
    freq
        If a requency calculation should be carried out.
    """

    name: str = "ORCA-Relax"
    orcasimpleinput: str = "Opt SlowConv NormalPrint"
    orcablock: str = ""
    freq: bool = False

    @job
    def make(self, atoms: Atoms) -> Dict[str, Any]:
        """
        Make the run.

        Parameters
        ----------
        atoms
            .Atoms object

        Returns
        -------
        Dict
            Summary of the run.
        """
        # TODO: Use " NumFreq" if analyticals aren't available.
        if (
            self.freq
            and "freq" not in self.orcasimpleinput.lower()
            and "numfreq" not in self.orcasimpleinput.lower()
        ):
            self.orcasimpleinput += " Freq"

        atoms.calc = ORCA(
            orcasimpleinput=self.orcasimpleinput,
            orcablock=self.orcablock,
        )
        atoms = run_calc(atoms)
        summary = summarize_run(atoms, ".out", additional_fields={"name": self.name})

        return summary

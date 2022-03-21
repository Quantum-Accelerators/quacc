import multiprocessing
from dataclasses import dataclass
from typing import Any, Dict

from ase.atoms import Atoms
from ase.calculators.orca import ORCA
from jobflow import Maker, job

from quacc.schemas.cclib import summarize_run
from quacc.util.basics import merge_dicts
from quacc.util.calc import run_calc


@dataclass
class StaticJob(Maker):
    """
    Class to carry out a single-point calculation.

    Parameters
    ----------
    name
        Name of the job.
    xc
        Exchange-correlation functional
    basis
        Basis set
    input_swaps
        Dictionary of orcasimpleinput swaps for the calculator.
        To enable new entries, set the value as True.
        To remove entries from the defaults, set the value as None/False.
    block_swaps
        Dictionary of orcablocks swaps for the calculator.
        To enable new entries, set the value as True.
        To remove entries from the defaults, set the value as None/False.
    """

    name: str = "ORCA-Static"
    xc: str = "wb97x-d3bj"
    basis: str = "def2-tzvp"
    input_swaps: Dict[str, Any] = None
    block_swaps: Dict[str, Any] = None

    @job
    def make(
        self, atoms: Atoms, charge: int = None, mult: int = None
    ) -> Dict[str, Any]:
        """
        Make the run.

        Parameters
        ----------
        atoms
            .Atoms object
        charge
            Charge of the system. If None, this is determined from the sum of
            atoms.get_initial_charges().
        mult
            Multiplicity of the system. If None, this is determined from 1+ the sum
            of atoms.get_initial_magnetic_moments().

        Returns
        -------
        Dict
            Summary of the run.
        """
        input_swaps = self.input_swaps or {}
        block_swaps = self.block_swaps or {}
        if not any(k for k in block_swaps if "nprocs" in k.lower()):
            nprocs = multiprocessing.cpu_count()
            block_swaps[f"%pal nprocs {nprocs} end"] = True

        default_inputs = {
            self.xc: True,
            self.basis: True,
            "sp": True,
            "slowconv": True,
            "normalprint": True,
        }
        default_blocks = {}

        inputs = merge_dicts(
            default_inputs, input_swaps, remove_none=True, remove_false=True
        )
        blocks = merge_dicts(
            default_blocks, block_swaps, remove_none=True, remove_false=True
        )
        orcasimpleinput = " ".join(list(inputs.keys()))
        orcablocks = " ".join(list(blocks.keys()))

        atoms.calc = ORCA(
            charge=charge if charge else round(sum(atoms.get_initial_charges())),
            mult=mult if mult else round(1 + sum(atoms.get_initial_magnetic_moments())),
            orcasimpleinput=orcasimpleinput,
            orcablocks=orcablocks,
        )
        atoms = run_calc(atoms)
        summary = summarize_run(
            atoms, "orca.out", additional_fields={"name": self.name}
        )

        return summary


@dataclass
class RelaxJob(Maker):
    """
    Class to carry out a geometry optimization.

    Parameters
    ----------
    name
        Name of the job.
    xc
        Exchange-correlation functional
    basis
        Basis set
    freq
        If a requency calculation should be carried out.
    input_swaps
        Dictionary of orcasimpleinput swaps for the calculator.
        To enable new entries, set the value as True.
        To remove entries from the defaults, set the value as None/False.
    block_swaps
        Dictionary of orcablocks swaps for the calculator.
        To enable new entries, set the value as True.
        To remove entries from the defaults, set the value as None/False.
    """

    name: str = "ORCA-Relax"
    xc: str = "wb97x-d3bj"
    basis: str = "def2-tzvp"
    freq: bool = False
    input_swaps: Dict[str, Any] = None
    block_swaps: Dict[str, Any] = None

    @job
    def make(
        self, atoms: Atoms, charge: int = None, mult: int = None
    ) -> Dict[str, Any]:
        """
        Make the run.

        Parameters
        ----------
        atoms
            .Atoms object
        charge
            Charge of the system. If None, this is determined from the sum of
            atoms.get_initial_charges().
        mult
            Multiplicity of the system. If None, this is determined from 1+ the sum
            of atoms.get_initial_magnetic_moments().

        Returns
        -------
        Dict
            Summary of the run.
        """
        input_swaps = self.input_swaps or {}
        block_swaps = self.block_swaps or {}
        if not any(k for k in block_swaps if "nprocs" in k.lower()):
            nprocs = multiprocessing.cpu_count()
            block_swaps[f"%pal nprocs {nprocs} end"] = True

        default_inputs = {
            self.xc: True,
            self.basis: True,
            "opt": True,
            "slowconv": True,
            "normalprint": True,
            "freq": True if self.freq else None,
        }
        default_blocks = {}

        inputs = merge_dicts(
            default_inputs, input_swaps, remove_none=True, remove_false=True
        )
        blocks = merge_dicts(
            default_blocks, block_swaps, remove_none=True, remove_false=True
        )
        orcasimpleinput = " ".join(list(inputs.keys()))
        orcablocks = " ".join(list(blocks.keys()))

        atoms.calc = ORCA(
            charge=charge if charge else round(sum(atoms.get_initial_charges())),
            mult=mult if mult else round(1 + sum(atoms.get_initial_magnetic_moments())),
            orcasimpleinput=orcasimpleinput,
            orcablocks=orcablocks,
        )
        atoms = run_calc(atoms)
        summary = summarize_run(
            atoms, "orca.out", additional_fields={"name": self.name}
        )

        return summary
